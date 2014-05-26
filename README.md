Enhancer Element Locator
========================

Tool for finding evolutionarily conserved mammalian enhancer
elements. For description of the underlying alignment algorithm, see
[Hallikas et.al. 2006](http://www.cs.helsinki.fi/u/kpalin/HallikasEtAl06.pdf) 

Some user guidance is available at
[http://www.cs.helsinki.fi/u/kpalin/EEL/] or in
[Palin et.al. 2006](https://www.cs.helsinki.fi/u/ukkonen/NatProtoc2006.pdf)


If you use EEL in academic publications, please cite:

> Hallikas, Palin,Sinjushina, Rautiainen, Partanen Ukkonen, Taipale: 
> Genome-wide Prediction of Mammalian Enhancers Based on Analysis of
> Transcription Factor Binding Affinity. 
> CELL 124, January 13, 2006.


Copyright Kimmo.Palin at helsinki.fi
Licence GPL

TF binding motif comparison by KL-divergence
--------------------------------------------

Newer versions of EEL contain a new command "computeKLdistances" which
computes all pairwise comparisons of loaded TF binding motifs. The
"help" command will explain the syntax. The method is used and
described in
[Wei et.al. 2010](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2905244/).
If you use this method in academic publications, please cite:

> Wei, G. H., Badis, G., Berger, M. F., Kivioja, T., Palin, K., Enge,
> M., ... & Taipale, J. (2010). Genome‐wide analysis of ETS‐family
> DNA‐binding in vitro and in vivo. The EMBO journal, 29(13),
> 2147-2160.



The article describes the method by:

> Comparison of binding profiles was performed using a novel algorithm
> that determines the similarity between TF motifs using the minimum
> Kullback–Leibler divergence between all translations and reverse
> complementations of the multinomial distributions defined by the
> motifs. Conceptually, the TF-motif divergence measures the information
> gained about the DNA sequence by knowledge of having binding sites for
> both of the two factors. The TF-motif divergence is defined as the
> minimum Kullback–Leibler divergence between all translations and
> reverse complementations of the multinomial distributions defined by
> the two TF motifs. The longer motif is inserted to a sequence with
> background distribution and the shorter motif is slid over the
> background/longer motif sequence. The KL divergence is computed
> between the multinomial distributions defined by (1) the shorter motif
> and (2) the part of the background/longer motif sequence overlapping
> the shorter motif. The same is repeated with the background/long motif
> sequence reverse complemented and the minimum of the KL divergences is
> taken. The TF-motif divergence is symmetric but does not fulfill the
> triangle inequality and thus is not a metric in the mathematical
> sense.
