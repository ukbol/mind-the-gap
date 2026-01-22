# [mind-the-gap](https://youtu.be/QExoX4ls9OM?si=EmPShkAIIfieSnqX)

## GenBank file batch downloader (edirect_genbank_fetch.sh)
Batch download all rbcL (ribulose-1,5-bisphosphate carboxylase) GenBank records from NCBI using EDirect. With a few changes to the script, it can be used to fetch genbank records for any target sequence.

**Requirements**
- EDirect (`conda install -c bioconda entrez-direct`)

**Usage**
```bash
sbatch edirect_genbank_fetch.sh
```
=

**Notes**
- 3-second delay between batches to respect NCBI rate limits
- Set `NCBI_API_KEY` environment variable for faster downloads (10 vs 3 requests/sec)
- Edit `QUERY` variable to modify search terms, or `BATCH` to change batch size.
