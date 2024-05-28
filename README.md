# Producing a table of DDX6 binding sites located within 3'UTRs across the transcriptome 
## Assessment specifications: 
_'You should take your own code or relevant software as a starting point, and modify or extend it so that it:_

_- is a github repository, with issues + pull requests that fully document the development process_

_- has unit testing with a good test coverage_

_- is pip installable_

_- is fully documented, e.g. with a README file and a "Read the Docs" website_

_- uses continuous integration (e.g. Github Actions) for automatic testing and documentation generation'_

## Background
### Argonaute Phosphorylation and the RNA-Induced Silencing Complex (RISC)
DDX6 is a RNA binding protein (RBP) that is recruited to the RNA-Induced Silencing Complex (RISC) following phosphorylation of the Argonaute (Ago2) protein. The RISC complex is responsible for seeking out specific messenger RNA transcripts (mRNAs) in the
cytoplasm and blocking their translation; Ago2 performs a key role in this function by guiding the RISC complex to target mRNA transcripts, or by directly cleaving mRNAs upon recognition. Ago2 is able to guide RISC to target mRNA transcripts with the help of
short RNA sequences called microRNAs (miRNAs). miRNA sequences are often complementary to 'seed sequences' located in the 3'-untranslated regions (3'UTRs) of certain mRNA transcripts. The loading of a single microRNA onto Ago2 prompts RISC assembly and a subsequent search for 
'seed sequences'. In this way, gene expression can be silenced. Recruitment of DDX6 to this complex following Ago2 phosphorylation is thought to help confer greater specificity for target mRNAs. This is because although the function of DDX6 
is relatively undefined, it is hypothesised to recognise specific secondary structures in RNA; therefore, it may contribute to target specificity in RISC function by binding to mRNA 3'UTRs with specific secondary structures.

### Accessing DDX6 binding data
We want to identify which mRNAs are likely to be targeted by RISC following Ago2 phosphorylation. One way in which this could be done is to identify targets for DDX6 binding. Fortunately, a series of enhanced crosslinking-immunoprecipitation experiments (eCLIP) experiments were carried
out to identify the binding sites for most RBPs across the trancsriptome as part of the ENCODE project [1]. Another database called POSTAR3 was produced using these publicly-available data, alongside data from similar experiments unaffiliated with the ENCODE project, to generate a 'comprehensive map' detailing the 
binding patterns of RBPs[2].


[1] https://pubmed.ncbi.nlm.nih.gov/32252787/
[2] https://academic.oup.com/nar/article/50/D1/D287/6353804?login=true
