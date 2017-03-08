#!/usr/bin/env nextflow
// Create fastq channel from samples.tsv

input_vcf = Channel.fromPath(params.vcf_file).map { file -> [file.baseName, false, file] }

input_vcf.into{input_vcf_plink; input_vcf_basicstats; input_vcf_popstats}

include_samples = params.vcfsamplenames
include_samples.removeAll(params.excludesamples)

sample_groups = []
file("${params.sample_groups}").readLines().each{line ->
  (sampleID, populationID) = line.split(/\t/)
  sample_groups.push([ "p":populationID, "s":sampleID ])
}

sample_groups = sample_groups.groupBy({it -> it.p}).collectEntries{[it.key, it.value.s]}
sample_groups = sample_groups.collect { key, value -> [key, value] }
sample_groups.each{a -> 
  a[1].removeAll(params.excludesamples)
  a[1].retainAll(params.vcfsamplenames)
}
sample_groups.removeAll({ it[1].empty })
sample_groups_file = file(params.sample_groups)

// This invalidates cache for some reason...
//
//
//
//
//

process Plink_bed {
  publishDir 'outputs/plink', mode: 'copy'
  tag {prefix}
  cpus 1
  memory 4.GB
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, pruned, file(vcf) from input_vcf_plink

  output:
  set prefix, pruned, file("*.bed"), file("*.bim"), file("*.fam") into plink_bed

  script:
  remove_file = file("temp-remove.txt")
  remove_file.text = "temp temp\n"
  params.excludesamples.each{
    remove_file.append("${it} ${it}\n")
  }
  """
  /usr/local/bin/plink --make-bed --remove ${remove_file} --vcf ${vcf} --allow-extra-chr --out temp
  /usr/local/bin/plink --make-bed --set-missing-var-ids @:#\\\$1,\\\$2 --allow-extra-chr --bfile temp --out ${prefix}
  rm temp.bed temp.bim temp.fam
  """
}

plink_bed.into{plink_bed_for_pruning; plink_bed}

process Plink_ld_pruning {
  publishDir 'outputs/plink', mode: 'copy'
  tag {prefix}
  cpus 1
  memory 4.GB
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, pruned, file(bed), file(bim), file(fam) from plink_bed_for_pruning

  output:
  set prefix, val(true), file("*.bed"), file("*.bim"), file("*.fam") into plink_pruned_bed
  set prefix, file("*.prune.in"), file("*.prune.out") into plink_pruned
  
  """
  /usr/local/bin/plink --indep 50 5 2 --allow-extra-chr --bed ${bed} --bim ${bim} --fam ${fam} --out ${prefix}
  /usr/local/bin/plink --make-bed --extract ${prefix}.prune.in --allow-extra-chr --bed ${bed} --bim ${bim} --fam ${fam} --out ${prefix}.ldpruned
  """
}

// Combine pruned and non-pruned genotypes
plink_bed = plink_bed.mix(plink_pruned_bed)
plink_bed.into{plink_bed_flat; plink_bed_stats}

plink_basicstats =
  [['param':'freq'       , 'ext':'frq'],
  ['param':'freq counts' , 'ext':'frq.counts'],
  ['param':'ibc'         , 'ext':'ibc'],
  ['param':'het'         , 'ext':'het'],
  ['param':'genome'      , 'ext':'genome'],
  ['param':'pca'         , 'ext':'eigen*'],
  ['param':'cluster'     , 'ext':'cluster*'],
  ['param':'homozyg'     , 'ext':'hom']]

process Plink_stats {
  publishDir 'outputs/plink_stats', mode: 'copy'
  tag {task.attempt + "." + plink_basicstat.param}
  cpus 1
  memory 4.GB
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, pruned, file(bed), file(bim), file(fam) from plink_bed_stats
  each plink_basicstat from plink_basicstats

  output:
  set prefix, pruned, val("${plink_basicstat.param}"), file("*.${plink_basicstat.ext}") into plink_basicstats_outputs

  """
  /usr/local/bin/plink --${plink_basicstat.param} --allow-extra-chr --bed ${bed} --bim ${bim} --fam ${fam} --out ${prefix}${pruned ? '-ldp' : '' }
  """
}

process Plink_traw {
  publishDir 'outputs/plink', mode: 'copy'
  tag {prefix + (pruned ? '-ldp' : '' )}
  cpus 1
  memory 4.GB
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, pruned, file(bed), file(bim), file(fam) from plink_bed_flat

  output:
  set prefix, pruned, file("*.traw") into plink_traw

  """
  /usr/local/bin/plink --recode A-transpose --allow-extra-chr --bed ${bed} --bim ${bim} --fam ${fam} --out ${prefix}${pruned ? '-ldp' : '' }
  """
}

plink_traw.into{plink_traw_idx; plink_traw_venn; plink_traw}

process Rf_idx {
  publishDir 'outputs/rarefaction', mode: 'copy'
  tag {prefix + (pruned ? '-ldp' : '' )}
  cpus 4
  memory 16.GB
  time 1.h
  errorStrategy { task.exitStatus == (143 | 139) ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set prefix, pruned, file(traw) from plink_traw_idx

  output:
  set prefix, pruned, file("*.Rdata") into rarefaction_index

  """
  rarefaction_index.R "SPACEHOLDER" 100 4 ${sample_groups_file} ${traw} ${prefix}${pruned ? '-ldp' : '' }
  """
}

//plink_traw = plink_traw.view{"plink_traw: $it"}
//rarefaction_index = rarefaction_index.view{"rf_idx: $it"}

rf_idx_traw = plink_traw.mix(rarefaction_index).groupTuple(by: [0,1]).map{a,b,c -> [a,b,c[0],c[1]]}
//rf_idx_traw = rf_idx_traw.view()
process Rf {
  publishDir 'outputs/rarefaction', mode: 'copy'
  tag {prefix + (pruned ? '-ldp' : '' )}
  cpus 8
  memory 32.GB
  time {task.attempt == 1 ? 6.h: 24.h}
  errorStrategy { task.exitStatus == (143 | 139) ? 'retry' : 'finish' }
  maxRetries 1
  maxErrors '-1'

  input:
  set prefix, pruned, file(traw), file(rfi) from rf_idx_traw
  each task_number from 1..params.rf_tasks
  val ntasks from params.rf_tasks

  output:
  set prefix, pruned, ntasks, file("*.Rdata") into rarefaction_batches

  """
  rarefaction.R ${task_number} ${ntasks} ${task.cpus} ${traw} ${rfi} ${prefix}${pruned ? '-ldp' : '' }
  """
}
//rarefaction_data = rarefaction_data.view()
rarefaction_batches = rarefaction_batches.groupTuple(by: [0,1,2])
//rarefaction_batches = rarefaction_batches.view()

process Rf_plots {
  publishDir 'outputs/rarefaction', mode: 'copy'
  tag {prefix + (pruned ? '-ldp' : '' )}
  cpus 1
  memory 16.GB
  time 2.d
  errorStrategy { task.exitStatus == (143 | 139) ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, pruned, ntasks, file(data) from rarefaction_batches

  output:
  set prefix, pruned, file("*.png"), file("*.pdf"), file("*.rf_concat.Rdata") into rarefaction_plots

  """
  rarefaction_plots.R ${ntasks} ${prefix}${pruned ? '-ldp' : '' }
  """
}


vcftools_basicstats =
   [['param':'counts'         , 'ext':'frq.count'],
    ['param':'depth'          , 'ext':'idepth'],
    ['param':'site-depth'     , 'ext':'ldepth'],
    ['param':'site-mean-depth', 'ext':'ldepth.mean'],
    ['param':'geno-depth'     , 'ext':'gdepth'],
    ['param':'het'            , 'ext':'het'],
    ['param':'relatedness'    , 'ext':'relatedness'],
    ['param':'relatedness2'   , 'ext':'relatedness2'],
    ['param':'missing-indv'   , 'ext':'imiss'],
    ['param':'singletons'     , 'ext':'singletons'],
    ['param':'freq'           , 'ext':'frq'],
    ['param':'site-quality'   , 'ext':'lqual']]

process BasicStats {
  publishDir 'outputs/basic_stats', mode: 'copy'
  tag {prefix + "-all-" + basicstat.param}

  cpus 1
  memory {task.attempt == 1 ? 4.GB: 8.GB}
  time {task.attempt == 1 ? 6.h: 24.h}
  errorStrategy { task.exitStatus == (143 | 139) ? 'retry' : 'finish' }
  maxRetries 1
  maxErrors '-1'

  input:
  set prefix, pruned, file(vcf) from input_vcf_basicstats
  each basicstat from vcftools_basicstats

  output:
  set prefix, pruned, val("${basicstat.param}"), file("*.${basicstat.ext}") into basic_stats_outputs

"""
/usr/local/bin/vcftools --vcf ${vcf} --${basicstat.param} --out ${prefix}${pruned ? '-ldp' : '' }-all
"""
}

basic_stats_outputs.into{basic_stats_outputs_frq; basic_stats_outputs_sq; basic_stats_outputs_smd}

vcftools_popstats =
   [['param':'counts'         , 'ext':'frq.count'],
    ['param':'depth'          , 'ext':'idepth'],
    ['param':'site-depth'     , 'ext':'ldepth'],
    ['param':'site-mean-depth', 'ext':'ldepth.mean'],
    ['param':'geno-depth'     , 'ext':'gdepth'],
    ['param':'het'            , 'ext':'het'],
    ['param':'freq'           , 'ext':'frq']]

process PopStats {
  publishDir 'outputs/pop_stats', mode: 'copy'
  tag {prefix + "-" + sample_group[0] + "-" + popstat.param}

  cpus 1
  memory {task.attempt == 1 ? 4.GB: 8.GB}
  time {task.attempt == 1 ? 6.h: 24.h}
  errorStrategy { task.exitStatus == (143 | 139) ? 'retry' : 'finish' }
  maxRetries 1
  maxErrors '-1'

  input:
  set prefix, pruned, file(vcf) from input_vcf_popstats
  each sample_group from sample_groups
  each popstat from vcftools_popstats

  output:
  set prefix, pruned, val("${sample_group[0]}"), val("${popstat.param}"), file("*.${popstat.ext}") into pop_stats_outputs

  script:
  keep_samples = sample_group[1].collect{"--indv $it"}.join(' ')

"""
/usr/local/bin/vcftools --vcf ${vcf} ${keep_samples} --${popstat.param} --out ${prefix}${pruned ? '-ldp' : '' }-${sample_group[0]}
"""
}

pop_stats_outputs.into{pop_stats_outputs_frq; pop_stats_outputs_sq; pop_stats_outputs_smd}
pop_stats_outputs_frq.into{pop_stats_outputs_frq_lgmaf; pop_stats_outputs_frq_maf}

//pop_stats_outputs_frq_lgmaf = pop_stats_outputs_frq_lgmaf.view()

pop_frq_channel_grouped = pop_stats_outputs_frq_lgmaf.filter{it[3] == "freq"}.groupTuple(by: [0,1]).map{a,b,c,d,e -> [a,b,c,e]}

process Venn {
  publishDir 'outputs/plots', mode: 'copy'

  cpus 8
  memory {task.attempt == 1 ? 32.GB: 64.GB}
  time {task.attempt == 1 ? 6.h: 24.h}
  errorStrategy { task.exitStatus == (143 | 139) ? 'retry' : 'finish' }
  maxRetries 1
  maxErrors '-1'

  input:
  set prefix, pruned, file(traw) from plink_traw_venn

  output:
  set prefix, pruned, file("*.pdf"), file("*.png") into venn_diagram_out

"""
plot_venn.R ${task.cpus} sumatrae,tigris,altaica,jacksoni ${sample_groups_file} ${traw} ${prefix}${pruned ? '-ldp' : '' }
"""
}

all_frq_channel = basic_stats_outputs_frq.filter{it[3] == "freq"}.map{a,b,c,d,e -> [a,b,c,e]}
pop_frq_channel = pop_stats_outputs_frq_maf.filter{it[3] == "freq"}.map{a,b,c,d,e -> [a,b,c,e]}

frq_channel = pop_frq_channel.mix(all_frq_channel)

process PlotMAF{
  publishDir 'outputs/plots', mode: 'copy'

  cpus 1
  memory {task.attempt == 1 ? 32.GB: 64.GB}
  time {task.attempt == 1 ? 6.h: 24.h}
  errorStrategy { task.exitStatus == (143 | 139) ? 'retry' : 'finish' }
  maxRetries 1
  maxErrors '-1'

  input:
  set prefix, pruned, population, file(frq) from frq_channel

  output:
  file("*.png") into plot_maf_png
  file("*.pdf") into plot_maf_pdf

  """
  plot_maf.R ${prefix}${pruned ? '-ldp' : '' }-${population} ${frq}
  """
}


/*all_sq_channel = basic_stats_outputs_sq.filter{it[2] == "site-quality"}.map{a,b,c,d -> [a,b,d]}
//all_sq_channel = all_sq_channel.view()
all_smd_channel = basic_stats_outputs_smd.filter{it[2] == "site-mean-depth"}.map{a,b,c,d -> [a,b,d]}
//all_smd_channel = all_smd_channel.view()
pq_channel = all_sq_channel.phase(all_smd_channel).view()
*/
/*
process PlotQUAL{
  publishDir 'outputs/plots', mode: 'copy'

  memory { 32.GB * task.attempt}
  time { 24.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 1
  maxErrors '-1'

  input:
  set prefix, pruned, file(sq), file(smd) from pq_channel

  output:
  file("*.png") into plot_qual_png
  file("*.pdf") into plot_qual_pdf

  """
  plot_qual_depth.R ${prefix}${pruned ? '-ldp' : '' } ${sq} ${smd}
  """
}
*/
/*process Vcfstats {
  publishDir 'outputs', mode: 'copy'
  tag {prefix}
  cpus 1
  memory 4.GB
  time 2.d
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set file(vcf), filterset from vcfstats_filtered_vcf

  output:
  set file("*.stats"), filterset into filtered_stats

  """
  /usr/local/opt/vcflib/bin/vcfstats ${vcf} > ${params.output_prefix}-snp-q${filterset.minqual}-gq${filterset.mingq}.stats
  """
}
*/


/*process Vcf2tsv {
  publishDir 'outputs', mode: 'copy'
  tag {prefix}
  cpus 1
  memory 4.GB
  time 2.d
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set file(vcf), filterset from vcf2tsv_filtered_vcf

  output:
  set file("*.tsv"), filterset into filtered_tsv

  """
  /usr/local/opt/vcflib/bin/vcf2tsv -n "NA" -g ${vcf} > ${params.output_prefix}-snp-q${filterset.minqual}-gq${filterset.mingq}.tsv
  """
}
*/

workflow.onComplete {
  println "Pipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
