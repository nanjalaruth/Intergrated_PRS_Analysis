nextflow.enable.dsl=2

process download_score_files {
    tag "Downloading score files from PGS catalogue"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(pheno)

    output:
        tuple val(blood_trait), val(pgs_id), path("${pgs_id}.txt"), path("${pgs_id}_metadata_scores.csv"), path(pheno)

    script:
        """
        curl -O https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pgs_id}/ScoringFiles/${pgs_id}.txt.gz
        gunzip ${pgs_id}.txt.gz
        
        curl -O https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pgs_id}/Metadata/${pgs_id}_metadata_scores.csv
        """
        
}

process conc_metadata {
    publishDir "${params.outdir}/${blood_trait}", mode: 'copy', overwrite: true
    tag "Concatenate metadata"

    input:
        tuple val(blood_trait), path(metadata_scores)

    output:
        tuple val(blood_trait), path("${blood_trait}_metadata_scores.csv")

    script:
        template "merge.R"

}

process modify_score_file {
    tag "Modifying score file from PGS catalogue"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(score_files), path(mdata), path(pheno)

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_edit.txt"), path(mdata), path(pheno)

    script:
        """
        #remove header information
        grep -v '^#' ${score_files} > ${blood_trait}_${pgs_id}_edit.txt
        """
}

process modify_score_file_2 {
    tag "Modifying score file from PGS catalogue"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(score_files), path(mdata), path(pheno), path(ref)

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_edit_2.txt"), path(mdata), path(pheno)

    script:
        
        """
        # Check if the columns chr_name and chr_position exist in the input file
        if awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++){if(\$i=="chr_name"){c1=1}if(\$i=="chr_position"){c2=1}}}END{exit !(c1 && c2)}'\
         ${score_files}; then \
             awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) { if (\$i == "chr_name") \
              { chr_name_index=i } else if (\$i == "chr_position") \
                { chr_position_index=i } } print \$0, "\\tSNP"; next} \
                 {print \$0 "\\t" \$chr_name_index ":" \$chr_position_index}' \
                   ${score_files} > ${blood_trait}_${pgs_id}_edit_2.txt
        else
            # Replace rsID with chr_name:chr_position using the reference file

            awk 'NR==FNR ? a[\$1] : \$1 in a' ${score_files} ${ref} | awk '{print \$1 "\\t" \$2":"\$3}' > ${blood_trait}_${pgs_id}_tmp
            join -t \$'\\t' ${score_files} ${blood_trait}_${pgs_id}_tmp | sed 's/chr_name:chr_position/SNP/g' > ${blood_trait}_${pgs_id}_edit_2.txt
        fi
        """
}

process modify_score_file_3 {
    tag "Modifying score file from PGS catalogue"
    
    input:
        tuple val(blood_trait), val(pgs_id), path(score_files), path(mdata), path(pheno)

    output:
        tuple val(blood_trait), val(pgs_id), path("${blood_trait}_${pgs_id}_modified.txt"), path(mdata),path(pheno)

    script:
        template "modify_score_file.R"
}


process compute_pgs_scores {
    tag "Computing PGS scores in PLINK"
    publishDir "${params.outdir}/${dataset}/${blood_trait}", mode: 'copy', overwrite: true
    
    input:
        tuple val(blood_trait), val(pgs_id), path(modified_scoring_file), path(mdata), path(pheno), val(dataset), path(bed), path(bim), path(fam)

    output:
        tuple val(dataset), val(blood_trait), val(pgs_id), path("${dataset}_${blood_trait}_${pgs_id}_prsval*")

    script:
        base = bed.baseName
        """
        /Users/rnanjala/Desktop/Polygenic_Risk_Scores/plink --bfile ${base} \\
        --score ${modified_scoring_file} 1 2 3 \\
        --pheno ${pheno} \\
        --out ${dataset}_${blood_trait}_${pgs_id}_prsval
        """
}


process modify_pgs_scores {
    tag "Modifying score file from PGS catalogue"
    
    input:
        tuple val(dataset), val(blood_trait), val(pgs_id), path(pgs_score)

    output:
        tuple val(dataset), val(blood_trait), val(pgs_id), path("${dataset}_${blood_trait}_${pgs_id}_pgs_scores.txt")

    script:
        """
         #extract FID and SCORE columns
        awk '{print \$1 "\\t" \$6}' ${pgs_score} > ${dataset}_${blood_trait}_${pgs_id}_pgs.txt
        awk -v OFS='\\t' '{if (NR==1) {\$2="${pgs_id}_SCORE"}; print}' ${dataset}_${blood_trait}_${pgs_id}_pgs.txt > ${dataset}_${blood_trait}_${pgs_id}_pgs_scores.txt
        """
        
}

process conc_scores {

    publishDir "${params.outdir}/${blood_trait}", mode: 'copy', overwrite: true
    
    tag "Concatenating pgs scores"
    
    input:
        tuple val(blood_trait), path(score)

    output:
        tuple val(blood_trait), path("${blood_trait}_pgs_scores.txt")

    script:
        template "merge_scores.R"
}

process intergrate_scores {
    tag "Modifying score file from PGS catalogue"
    publishDir "${params.outdir}/${blood_trait}", mode: 'copy', overwrite: true
    
    input:
        tuple val(blood_trait), path(score_files)
        tuple val(blood_trait), path(pheno), path(cov)

    output:
        tuple val(blood_trait), path("${blood_trait}_protsignature"),
         path("${blood_trait}_correlation.pdf"), path("${blood_trait}_pred_score.txt"),
         path("${blood_trait}_ug_prs_pcs.txt")

    script:
        template "intergrate_scores.R"
}

process association_analysis {
    tag "Association analysis"
    publishDir "${params.outdir}/${blood_trait}", mode: 'copy', overwrite: true
    
    input:
        tuple val(blood_trait), path(ug_prs_pcs), path(lipid_trait)

    output:
        tuple val(blood_trait), path("${blood_trait}_forest_plot.pdf")

    script:
        template "prediction.R"
}

workflow{
     //step 1
     // Download score files from PGS catalogue
    input = Channel.fromList(params.pheno)
    //input.view()
    download_score_files(input)

    //step 1.2
    inpt = download_score_files.out
        .map{ btrait, pgsid, score_file, mdata, pheno -> [btrait, mdata] }
        .groupTuple()
    //inpt.view()
    conc_metadata(inpt)

    //step 2
    // Modify score file
    modify_ch = download_score_files.out
        //.map{ btrait, pgsid, score_file, mdata, pheno -> 
          //          return [btrait, pgsid, score_file, pheno] 
            //    }
    //modify_ch.view()
    modify_score_file(modify_ch)

    //step 2.2
    ref = Channel.fromPath(params.reference)
    input = modify_score_file.out.combine(ref)
    //input.view()
    modify_score_file_2(input)

    //step 2.3
    modify_3_ch = modify_score_file_2.out
    //modify_3_ch.view()
    modify_score_file_3(modify_3_ch)

    //step 3
    // calculate pgs score
    modify_score_out = modify_score_file_3.out
    plink_ch = Channel.fromList(params.plink_file)
    pgs_ch = modify_score_out.combine(plink_ch)
    //pgs_ch.view()
    compute_pgs_scores(pgs_ch)

    //step 3.2 modify_pgs_scores
    in = compute_pgs_scores.out
     .map{dset, btrait, pgsid, pgs_score -> [dset, btrait, pgsid, pgs_score[3]]}
     //in.view()
    modify_pgs_scores(in)

    //step 3.3 conc scores
    input = modify_pgs_scores.out
        .map{dset, btrait, pgsid, pgs_score -> [btrait, pgs_score]}
        .groupTuple()
    //input.view()
    conc_scores(input)

    //step 4 Intergrate scores
    scores = conc_scores.out
    //scores.view()
    pheno = Channel.fromList(params.phenotype)
    cov = Channel.fromPath(params.covariate)
    cov_pheno = pheno
        .combine(cov)    
    //cov_pheno.view()
    intergrate_scores(scores, cov_pheno)

    //step 5 Predictions
    trait = Channel.fromPath(params.trait)
    input = intergrate_scores.out
        .map{btrait, sig, cor, pred, ug_prs_pcs -> [btrait, ug_prs_pcs]}
        .combine(trait)
    //input.view()
    association_analysis(input)
   
}