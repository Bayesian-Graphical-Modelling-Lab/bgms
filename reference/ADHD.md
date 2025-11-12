# ADHD Symptom Checklist for Children Aged 6–8 Years

This dataset includes ADHD symptom ratings for 355 children aged 6 to 8
years from the Children’s Attention Project (CAP) cohort (Silk et al.
2019) . The sample consists of 146 children diagnosed with ADHD and 209
without a diagnosis. Symptoms were assessed through structured
interviews with parents using the NIMH Diagnostic Interview Schedule for
Children IV (DISC-IV) (Shaffer et al. 2000) . The checklist includes 18
items: 9 Inattentive (I) and 9 Hyperactive/Impulsive (HI). Each item is
binary (1 = present, 0 = absent).

## Usage

``` r
data("ADHD")
```

## Format

A matrix with 355 rows and 19 columns.

- group:

  ADHD diagnosis: 1 = diagnosed, 0 = not diagnosed

- avoid:

  Often avoids, dislikes, or is reluctant to engage in tasks that
  require sustained mental effort (I)

- closeatt:

  Often fails to give close attention to details or makes careless
  mistakes in schoolwork, work, or other activities (I)

- distract:

  Is often easily distracted by extraneous stimuli (I)

- forget:

  Is often forgetful in daily activities (I)

- instruct:

  Often does not follow through on instructions and fails to finish
  schoolwork, chores, or duties in the workplace (I)

- listen:

  Often does not seem to listen when spoken to directly (I)

- loses:

  Often loses things necessary for tasks or activities (I)

- org:

  Often has difficulty organizing tasks and activities (I)

- susatt:

  Often has difficulty sustaining attention in tasks or play activities
  (I)

- blurts:

  Often blurts out answers before questions have been completed (HI)

- fidget:

  Often fidgets with hands or feet or squirms in seat (HI)

- interrupt:

  Often interrupts or intrudes on others (HI)

- motor:

  Is often "on the go" or often acts as if "driven by a motor" (HI)

- quiet:

  Often has difficulty playing or engaging in leisure activities quietly
  (HI)

- runs:

  Often runs about or climbs excessively in situations in which it is
  inappropriate (HI)

- seat:

  Often leaves seat in classroom or in other situations in which
  remaining seated is expected (HI)

- talks:

  Often talks excessively (HI)

- turn:

  Often has difficulty awaiting turn (HI)

## Source

Silk et al. (2019) . Data retrieved from
[doi:10.1371/journal.pone.0211053.s004](https://doi.org/10.1371/journal.pone.0211053.s004)
. Licensed under the CC-BY 4.0:
https://creativecommons.org/licenses/by/4.0/

## References

Shaffer D, Fisher P, Lucas CP, Dulcan MK, Schwab-Stone ME (2000). “NIMH
Diagnostic Interview Schedule for Children Version IV (NIMH DISC-IV):
description, differences from previous versions, and reliability of some
common diagnoses.” *Journal of the American Academy of Child &
Adolescent Psychiatry*, **39**, 28–38.
[doi:10.1097/00004583-200001000-00014](https://doi.org/10.1097/00004583-200001000-00014)
, PMID: 10638065.  
  
Silk TJ, Malpas CB, Beare R, Efron D, Anderson V, Hazell P, Jongeling B,
Nicholson JM, Sciberras E (2019). “A network analysis approach to ADHD
symptoms: More than the sum of its parts.” *PLOS ONE*, **14**(1),
e0211053.
[doi:10.1371/journal.pone.0211053](https://doi.org/10.1371/journal.pone.0211053)
.
