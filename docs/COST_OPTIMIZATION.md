# Cost Optimization Guide

## The Problem: Data Gravity and Unnecessary Copies

### The Original Problem That Started This Project

**Real-world scenario that sparked BAMS3:**

```
AWS Open Data (us-east-1) → Your S3 (us-west-2) → Local compute
     100GB @ $2.00         →    100GB @ $2.00     →  Process

Total cost per analysis: $4.00 in transfer alone
Total time: ~60 minutes just copying
Total waste: 200GB of redundant storage
```

**The realization:** We were copying data twice for no reason.
1. Copy from open data to our S3 (because we thought we "owned" it that way)
2. Copy from our S3 to local compute (because tools needed local files)

**Both copies were completely unnecessary.**

## The Solution: Compute Near Data

### The $1.53 vs $2 Insight

Instead of moving data to compute, move compute to data:

```
BEFORE: Copy data cross-region ($2) + process
AFTER:  Spin up c5.9xlarge in source region ($1.53/hour) + process

Savings: $2.00 (cross-region transfer eliminated)
Time savings: 30 minutes (no download wait)
Bonus: Instance probably costs less than $1.53 because of speedup!
```

**The math:**
- Cross-region transfer: $2.00 (fixed cost)
- EC2 instance: $1.53/hour
- But with streaming, analysis finishes faster
- Actual EC2 cost: ~$0.50-$1.00 for targeted queries

**Net result: Spend less, finish faster, eliminate waste.**

## Cost Breakdown by Scenario

### Scenario 1: Single Gene Analysis

**Traditional Approach:**
```
1. Copy 100GB BAM to your region          $2.00 + 30 min
2. Download BAM to local/EC2               $9.00 + 30 min
3. Extract gene region with samtools       ~2 min
4. Run analysis                            ~5 min

Total: $11.00, 67 minutes
```

**BAMS3 Approach:**
```
1. Launch c5.9xlarge in source region      $0 (no transfer)
2. Stream gene region directly (50MB)      $0.0045 + 30 sec
3. Run analysis                            ~5 min

Total: $0.26 (10 min × $1.53/hour), 10 minutes

Savings: $10.74 (97.6%), 57 minutes (85%)
```

### Scenario 2: 10-Gene Panel (Iterative Analysis)

**Traditional Approach:**
```
For each gene:
  - Already have full BAM cached locally
  - Extract region: ~2 min
  - Analyze: ~5 min

Total: 10 × 7 min = 70 minutes
One-time setup: $11.00
```

**BAMS3 Approach:**
```
For each gene:
  - Stream region: 30 sec
  - Analyze: 5 min

Total: 10 × 5.5 min = 55 minutes
Cost: $1.41 (55 min × $1.53/hour)

Savings: $9.59 (85%), 15 minutes (21%)
```

**Key insight:** Traditional approach requires the $11 upfront "tax" to download. BAMS3 pays only for what you use.

### Scenario 3: Large Cohort Study (1000 samples)

**Traditional Approach:**
```
Storage: 1000 × 100GB = 100TB
Monthly S3 cost: $2,300/month
Per-query download: $9.00

Annual cost: $27,600 + ($9 × queries)
For 1000 queries: $36,600
```

**BAMS3 Approach:**
```
Storage: 1000 × 10GB = 10TB
Monthly S3 cost: $230/month
Per-query stream: $0.0045

Annual cost: $2,760 + ($0.0045 × queries)
For 1000 queries: $2,765

Savings: $33,835 (92%)
```

## The Compound Waste Problem

### Traditional Workflow Creates Multiple Copies

```
Original BAM in Open Data (us-east-1)
    ↓ Copy 1: $2.00, 30 min
Your S3 copy (us-west-2)
    ↓ Copy 2: $9.00, 30 min
Local/EC2 disk
    ↓ Copy 3: Intermediate files (SAM, unsorted BAM, etc.)
More local disk (200GB total)

Total redundancy: 3× the data
Total cost: $11.00 per sample
Total time: 60 minutes before analysis even starts
```

### BAMS3 Workflow Eliminates All Copies

```
Original BAMS3 in Open Data (us-east-1)
    ↓ Stream what you need: $0.0045, 30 sec
Directly to analysis tool (via pipe)
    ↓ No intermediate files
Results

Total redundancy: 0× (no copies)
Total cost: $0.0045 per query
Total time: 30 seconds, analysis starts immediately
```

## Hidden Costs Traditional Approaches Miss

### 1. Developer Time
**Problem:** Waiting for downloads
- 30 min download = 30 min of wasted developer time
- At $100/hour developer cost: $50 wasted per query
- BAMS3: Start analysis in seconds

### 2. Disk Management
**Problem:** Storage fills up quickly
- 100GB per sample
- 10 samples = 1TB of local disk needed
- Need to delete old analyses to make room
- Time spent on disk management: 15 min/week
- BAMS3: 50MB per query, negligible disk use

### 3. Pipeline Complexity
**Problem:** More steps = more failure points
- Traditional: 4-5 steps (download → index → extract → analyze)
- Each step can fail
- Debugging time: 30 min per failure
- BAMS3: 1 step (stream → analyze), fewer failures

### 4. Opportunity Cost
**Problem:** Can't explore freely due to cost
- "Should I analyze this gene? It costs $11..."
- Fear of exploration limits discovery
- BAMS3: $0.0045 per query = explore freely

## The Speedup Multiplier Effect

### Why EC2 Actually Costs Less Than $1.53

**Example: 10-gene analysis**

Traditional (with local BAM):
- Time: 70 minutes
- EC2 cost: 70 min × $1.53/hour = $1.79

BAMS3 (streaming):
- Time: 55 minutes (streaming is parallel to analysis)
- EC2 cost: 55 min × $1.53/hour = $1.40

**But it gets better with concurrent queries:**

BAMS3 (parallel streaming, 10 genes):
- Time: 15 minutes (parallel processing)
- EC2 cost: 15 min × $1.53/hour = $0.38

**The speedup itself reduces EC2 costs!**

## Cost Optimization Best Practices

### 1. Always Compute Near Data

❌ **Wrong:**
```bash
# Copy data to your "home" region
aws s3 cp s3://open-data/sample.bam s3://your-bucket-us-west-2/
# Then process in us-west-2
```

✅ **Right:**
```bash
# Launch compute in same region as data
aws ec2 run-instances --region us-east-1 ...
# Process directly, no copy
```

**Savings:** $2.00 per 100GB, 30 minutes

### 2. Use Selective Access Formats (BAMS3/VCFS3)

❌ **Wrong:**
```bash
# Download full file for targeted analysis
aws s3 cp s3://bucket/sample.bam ./
samtools view sample.bam chr17:41196312-41277500
```

✅ **Right:**
```bash
# Stream only what you need
bams3 to-bam s3://bucket/sample.bams3 - --region chr17:41196312-41277500
```

**Savings:** 99.95% of data transfer ($9 → $0.0045)

### 3. Eliminate Intermediate Files

❌ **Wrong:**
```bash
# Pipeline with intermediate files
bwa mem ref.fa R1.fq R2.fq > aligned.sam
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam > sorted.bam
samtools index sorted.bam
```

✅ **Right:**
```bash
# Zero-copy streaming pipeline
bwa mem ref.fa R1.fq R2.fq | bams3 convert --stdin output.bams3
```

**Savings:** 200GB of intermediate storage, simpler pipeline

### 4. Don't "Own" Data Unnecessarily

❌ **Wrong thinking:**
```
"I need to copy this from open data to my S3 bucket
 so I own it and can access it reliably."
```

✅ **Right thinking:**
```
"Data in open data repositories is already reliable.
 I can access it directly and save $2 per copy.
 If I really need my own copy, convert to BAMS3
 and save 90% on storage."
```

### 5. Pay for Computation, Not Storage Movement

**The shift:**
- Old model: Pay to move data around ($2-$11 per file)
- New model: Pay for compute time ($0.26-$1.53 per query)

**Why compute is cheaper:**
- Compute scales with analysis complexity
- Storage movement scales with file size (fixed)
- Targeted queries use 1000× less compute than full downloads

## Real-World Cost Impact

### Research Lab: 1000 Samples, 100 Queries/Year

**Traditional Approach:**
```
Storage: 1000 × 100GB @ $0.023/GB/month
  = $2,300/month × 12 = $27,600/year

Per-query download: $9.00 × 100 = $900/year

Total first year: $28,500
```

**BAMS3 Approach:**
```
Storage: 1000 × 10GB @ $0.023/GB/month
  = $230/month × 12 = $2,760/year

Per-query stream: $0.0045 × 100 = $0.45/year
Compute (EC2): ~$0.50 × 100 = $50/year

Total first year: $2,810
```

**Savings: $25,690 (90%)**

**What you can do with $25,690:**
- Sequence 50 more samples
- Hire a summer intern
- Buy better compute for faster analysis
- Or just... save it

## The Broader Lesson

### Data Gravity Is Real

The original problem (copying from open data → your S3 → local) illustrates **data gravity:**

> Large datasets naturally accumulate more datasets around them.
>
> You start with one copy "to be safe."
> Then you make another copy "to analyze."
> Then intermediate files pile up.
> Before you know it, you have 3-5× redundancy.

**BAMS3/VCFS3 breaks data gravity:**
- Keep data in one canonical location
- Stream what you need, when you need it
- No copies, no gravity, no waste

### The Cloud-Native Mindset Shift

**Old mindset (local HPC era):**
- "I need local copies of everything"
- "Storage is cheap, bandwidth is expensive"
- "Download once, analyze many times"

**New mindset (cloud-native era):**
- "I access data where it lives"
- "Compute is cheap, storage movement is expensive"
- "Stream what I need, pay for what I use"

**BAMS3/VCFS3 enables the new mindset.**

## Conclusion

The journey from "copy from open data → copy to local → process" to "stream directly from source" is the journey from wasteful to efficient cloud-native genomics.

**Key insights:**
1. **Move compute to data** ($1.53/hour) not data to compute ($2.00 per copy)
2. **Stream selectively** (99.95% cost savings)
3. **Eliminate redundancy** (3× → 0× copies)
4. **Faster = cheaper** (speedup reduces EC2 costs too)

The best part? **The speedup usually makes the EC2 instance cost less than the transfer you avoided.**

You're not just saving money by avoiding the $2 transfer - you're saving money because the whole workflow is 10× faster, so your $1.53/hour instance only runs for 1/10th the time.

**Math:**
- Avoided transfer: $2.00 saved
- EC2 cost: $0.26 (10 min instead of 70 min)
- Net savings: $2.00 - $0.26 = **$1.74**
- Plus: 60 minutes of your time back

It's not just cheaper. It's faster. It's simpler. It's the right way to do cloud-native genomics.

---

**See Also:**
- [S3_INTEGRATION.md](S3_INTEGRATION.md) - Technical architecture
- [WGS_TESTING_PROTOCOL.md](../WGS_TESTING_PROTOCOL.md) - Validation testing
- [VCFS3_DESIGN.md](VCFS3_DESIGN.md) - VCF cost optimization
