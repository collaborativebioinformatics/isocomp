

There seem to be three important bits of code surrounding fusion gene finding in 
(cupcake/tofu)[https://github.com/Magdoll/cDNA_Cupcake/tree/master/cupcake/tofu]

First, (fusion gene candidates are identified)[https://github.com/Magdoll/cDNA_Cupcake/blob/81b7e7f6aeb53e15c11dd30a68a498a58e5f390a/cupcake/tofu/fusion_finder.py#L225].

```python
def find_fusion_candidates(
  sam_filename, 
  query_len_dict, 
  min_locus_coverage=.05, 
  min_locus_coverage_bp=1, 
  min_total_coverage=.99, 
  min_dist_between_loci=10000, 
  min_identity=0.95):
   """
    Return dict of
       fusion candidate qID --> list (in order) of the fusion ranges (ex: (chr3,100,200), (chr1,500,1000))
    (1) must map to 2 or more loci
    (2) minimum coverage for each loci is 5% AND minimum coverage in bp is >= 1 bp
    (3) total coverage is >= 95%
    (4) distance between the loci is at least 10kb
    """
```

 Those candidates are filtered and subsequently reported 
[here](https://github.com/Magdoll/cDNA_Cupcake/blob/81b7e7f6aeb53e15c11dd30a68a498a58e5f390a/cupcake/tofu/fusion_finder.py#L262). 

```python

def fusion_main(
  fa_or_fq_filename, 
  sam_filename, 
  output_prefix, 
  cluster_report_csv=None,
  is_fq=False, 
  allow_extra_5_exons=True, 
  skip_5_exon_alt=True,
  min_locus_coverage=.05, 
  min_total_coverage=.99,
  min_locus_coverage_bp=1, 
  min_dist_between_loci=10000,
  min_identity=0.95,
  is_flnc=False):
 """
    (1) identify fusion candidates (based on mapping, total coverage, identity, etc)
    (2) group/merge the fusion exons, using an index to point to each individual part
    (3) use BranchSimple to write out a tmp GFF where
         PBfusion.1.1 is the first part of a fusion gene
         PBfusion.1.2 is the second part of a fusion gene
    (4) read the tmp file from <3> and modify it so that
         PBfusion.1 just represents the fusion gene (a single transcript GFF format)

fusion candidates are defined like so:
```

If there are multiple samples, both fusion and other isoforms may be chained 
together, which I understand as a `merge` step. Why this is called 'chaining' vs 
'merging', I'm not sure. 'merge' would be far clear-er, so maybe there is a subtlety 
that I am missing?

To 'chain' (merge) the non fusion genes, a [data structure built around an 
IntervalTree](https://github.com/Magdoll/cDNA_Cupcake/blob/81b7e7f6aeb53e15c11dd30a68a498a58e5f390a/cupcake/tofu/counting/combine_abundance_across_samples.py#L98) is used. The fusion chaining uses a 
[child](https://github.com/Magdoll/cDNA_Cupcake/blob/81b7e7f6aeb53e15c11dd30a68a498a58e5f390a/cupcake/tofu/counting/combine_abundance_across_samples.py#L280) of this structure. The goal is to 
use these structures to create a non redundant set of isoforms. In both cases, the 
constructors (links above are to the class definition) 
expose some settings which tune what is considered a match.