# RVIDS3: Research Video Format for Cloud Storage

## Overview

RVIDS3 (Research Video for S3) is a cloud-native format for behavioral research video with tightly integrated annotations, transcriptions, and temporal analysis. Designed for psychology, linguistics, education research, and behavioral sciences.

## Problem Statement

### Current Research Video Workflow

Researchers studying human behavior face significant challenges:

1. **Fragmented Data**
   - Video files: .mp4, .mov (large, monolithic)
   - Audio transcripts: .txt, .docx (separate files)
   - Annotations: ELAN .eaf, Anvil .anvil (proprietary XML)
   - Metadata: Excel spreadsheets (disconnected)
   - Analysis: Custom scripts (not reproducible)

2. **Storage and Access Issues**
   - Download entire 2-hour video to analyze 30 seconds of behavior
   - No selective access by annotation or speaker
   - Expensive to store and transfer
   - Difficult to share with collaborators

3. **Analysis Challenges**
   - "Find all instances of [behavior]" requires manual search
   - No programmatic access to time-aligned data
   - Difficult to aggregate across multiple sessions
   - Limited tool interoperability

### Example Use Cases

**Psychology Research**:
- Study parent-child interactions during play
- Extract all moments where child exhibits target behavior
- Analyze temporal patterns across 100+ hour-long sessions
- Share specific clips with collaborators

**Linguistics**:
- Analyze conversation dynamics in dyads
- Extract turns by speaker
- Study prosody, pauses, overlaps
- Time-align transcript with video/audio

**Education Research**:
- Classroom observation studies
- Extract teacher questions and student responses
- Analyze instructional sequences
- Code behaviors across multiple observers

## Design Goals

1. **Time-Based Selective Access**: Query by timestamp, duration, or time range
2. **Annotation-Based Query**: "Find all moments tagged as [behavior]"
3. **Speaker/Actor Query**: Extract all video segments for specific speaker
4. **Multi-Layer Annotations**: Support multiple annotation tiers (behaviors, speech, gestures)
5. **Cloud-Native**: Store and query directly from S3
6. **Tool Compatibility**: Export to ELAN, Anvil, VCode, standard video formats
7. **Cost Efficiency**: 90%+ savings on selective queries

## Architecture

### Format Structure

```
session.rvids3/
├── _metadata.json              # Dataset metadata
├── _manifest.json              # Segment and chunk manifest
├── _header.json                # Schema, tiers, speakers
├── _index/
│   ├── time.idx                # Temporal index (binary)
│   ├── annotation.idx          # Annotation tag index
│   ├── speaker.idx             # Speaker/actor index
│   └── tier.idx                # Annotation tier index
├── video/
│   ├── segments/
│   │   ├── 00000-00030.mp4     # 30-second video segments
│   │   ├── 00030-00060.mp4
│   │   └── ...
├── audio/
│   ├── segments/
│   │   ├── 00000-00030.opus    # 30-second audio segments (Opus codec)
│   │   ├── 00030-00060.opus
│   │   └── ...
├── transcript/
│   ├── words.json.zst          # Word-level timestamps
│   └── utterances.json.zst     # Utterance boundaries
└── annotations/
    ├── tiers/
    │   ├── behavior.json.zst   # Behavior annotation tier
    │   ├── gesture.json.zst    # Gesture annotation tier
    │   └── interaction.json.zst
    └── timeline.json.zst       # Master timeline

```

### Metadata Format

#### `_metadata.json`

```json
{
  "format": "RVIDS3",
  "version": "0.1.0",
  "session_id": "parent_child_001",
  "created": "2026-01-16T10:30:00Z",
  "duration_seconds": 3600,
  "video": {
    "width": 1920,
    "height": 1080,
    "fps": 30,
    "codec": "h264",
    "total_size_bytes": 2147483648
  },
  "audio": {
    "sample_rate": 48000,
    "channels": 2,
    "codec": "opus",
    "total_size_bytes": 67108864
  },
  "participants": [
    {
      "id": "parent_01",
      "role": "parent",
      "metadata": {"age": 35, "gender": "F"}
    },
    {
      "id": "child_01",
      "role": "child",
      "metadata": {"age": 4, "gender": "M"}
    }
  ],
  "annotation_tiers": [
    {
      "name": "behavior",
      "type": "point",
      "vocabulary": ["reach", "grasp", "look", "point", "vocalize"]
    },
    {
      "name": "interaction",
      "type": "interval",
      "vocabulary": ["joint_attention", "parallel_play", "conflict"]
    }
  ],
  "segmentation": {
    "video_segment_seconds": 30,
    "audio_segment_seconds": 30,
    "total_segments": 120
  }
}
```

#### `_manifest.json`

```json
{
  "segments": [
    {
      "index": 0,
      "start_time": 0.0,
      "end_time": 30.0,
      "video": {
        "path": "video/segments/00000-00030.mp4",
        "size_bytes": 18874368,
        "keyframes": [0.0, 1.0, 2.0, ...]
      },
      "audio": {
        "path": "audio/segments/00000-00030.opus",
        "size_bytes": 524288
      },
      "annotations": [
        {"tier": "behavior", "count": 5},
        {"tier": "interaction", "count": 2}
      ],
      "transcript": {
        "utterances": 3,
        "words": 42
      }
    },
    ...
  ]
}
```

### Annotation Format

#### `annotations/tiers/behavior.json.zst`

```json
{
  "tier": "behavior",
  "type": "point",
  "annotations": [
    {
      "id": "beh_001",
      "time": 5.234,
      "code": "reach",
      "participant": "child_01",
      "metadata": {
        "hand": "right",
        "target": "toy_car",
        "success": true
      }
    },
    {
      "id": "beh_002",
      "time": 5.891,
      "code": "grasp",
      "participant": "child_01",
      "metadata": {
        "hand": "right",
        "object": "toy_car"
      }
    },
    ...
  ]
}
```

#### `annotations/tiers/interaction.json.zst`

```json
{
  "tier": "interaction",
  "type": "interval",
  "annotations": [
    {
      "id": "int_001",
      "start_time": 10.5,
      "end_time": 45.2,
      "code": "joint_attention",
      "participants": ["parent_01", "child_01"],
      "metadata": {
        "initiator": "parent_01",
        "target_object": "picture_book"
      }
    },
    ...
  ]
}
```

### Transcript Format

#### `transcript/utterances.json.zst`

```json
{
  "utterances": [
    {
      "id": "utt_001",
      "speaker": "parent_01",
      "start_time": 3.45,
      "end_time": 5.67,
      "text": "Do you see the car?",
      "confidence": 0.96
    },
    {
      "id": "utt_002",
      "speaker": "child_01",
      "start_time": 6.12,
      "end_time": 6.89,
      "text": "Car!",
      "confidence": 0.91
    },
    ...
  ]
}
```

#### `transcript/words.json.zst`

```json
{
  "words": [
    {
      "utterance_id": "utt_001",
      "word": "Do",
      "start_time": 3.45,
      "end_time": 3.52,
      "confidence": 0.98
    },
    {
      "utterance_id": "utt_001",
      "word": "you",
      "start_time": 3.53,
      "end_time": 3.61,
      "confidence": 0.97
    },
    ...
  ]
}
```

### Index Formats

#### Time Index (`_index/time.idx`)

Binary index for fast temporal queries:

```
[Header: 32 bytes]
  - Magic: "RVID_TIME_IDX" (14 bytes)
  - Version: uint16 (2 bytes)
  - Num segments: uint32 (4 bytes)
  - Total duration: float64 (8 bytes)
  - Reserved: 4 bytes

[Segment entries: 64 bytes each]
  - Start time: float64 (8 bytes)
  - End time: float64 (8 bytes)
  - Video offset: uint64 (8 bytes)
  - Video size: uint64 (8 bytes)
  - Audio offset: uint64 (8 bytes)
  - Audio size: uint64 (8 bytes)
  - Annotation count: uint32 (4 bytes)
  - Transcript count: uint32 (4 bytes)
  - Reserved: 16 bytes
```

#### Annotation Index (`_index/annotation.idx`)

```json
{
  "tier": "behavior",
  "index": {
    "reach": [
      {"time": 5.234, "segment": 0},
      {"time": 67.89, "segment": 2},
      ...
    ],
    "grasp": [
      {"time": 5.891, "segment": 0},
      ...
    ]
  }
}
```

#### Speaker Index (`_index/speaker.idx`)

```json
{
  "speakers": {
    "parent_01": {
      "total_utterances": 234,
      "total_duration": 1023.45,
      "segments": [0, 1, 2, 3, 5, 7, ...]
    },
    "child_01": {
      "total_utterances": 156,
      "total_duration": 456.78,
      "segments": [0, 1, 2, 4, 6, 8, ...]
    }
  }
}
```

## Query Patterns

### 1. Time-Based Queries

**Extract 10-second clip starting at 1:30**:
```bash
rvids3 extract session.rvids3 output.mp4 \
  --start 90.0 --duration 10.0
```

**S3 Transfer**: 2 segments × 20MB = 40MB (vs 2GB full video)
**Savings**: 98%

### 2. Annotation-Based Queries

**Find all "joint_attention" moments**:
```bash
rvids3 query session.rvids3 \
  --tier interaction \
  --code joint_attention \
  --output clips/
```

**Result**: 15 clips (3-60 seconds each)
**S3 Transfer**: 50MB (vs 2GB full video)
**Savings**: 97.5%

### 3. Speaker/Participant Queries

**Extract all child vocalizations**:
```bash
rvids3 query session.rvids3 \
  --speaker child_01 \
  --tier behavior \
  --code vocalize \
  --context 2.0
```

**Result**: 42 clips with 2-second context
**S3 Transfer**: 30MB
**Savings**: 98.5%

### 4. Multi-Condition Queries

**Complex query: Parent questions followed by child response**:
```bash
rvids3 query session.rvids3 \
  --sequence \
    "parent_01:question" \
    "child_01:response" \
  --max-gap 5.0 \
  --output analysis/turn_taking/
```

### 5. Export to Standard Formats

**Export to ELAN**:
```bash
rvids3 export session.rvids3 output.eaf --format elan
```

**Export to Anvil**:
```bash
rvids3 export session.rvids3 output.anvil --format anvil
```

**Export specific tier to CSV**:
```bash
rvids3 export session.rvids3 behaviors.csv \
  --format csv \
  --tier behavior
```

## Compression and Segmentation

### Video Segmentation Strategy

**Segment Size**: 30 seconds (configurable)

**Why 30 seconds:**
- Typical behavioral units: 10-60 seconds
- GOP (Group of Pictures) alignment: 30s = 900 frames at 30fps
- Manageable download size: 15-25MB per segment
- Balance between granularity and overhead

**Codec**: H.264 with frequent keyframes
- Keyframe every 1 second (30 frames)
- Allows precise seeking within segments
- Standard codec with universal compatibility

### Audio Segmentation

**Segment Size**: 30 seconds (aligned with video)

**Codec**: Opus
- Superior quality at low bitrates (32-64 kbps)
- Low latency encoding
- Open standard
- ~500KB per 30-second segment

### Annotation Compression

**Format**: JSON with zstd compression
- Human-readable base format
- High compression ratio (10:1 typical)
- Fast decompression
- Streaming support

**Typical Sizes** (1-hour session):
- Raw annotations: 2MB
- Compressed: 200KB

## Implementation Plan

### Phase 1: Core Format and Converter (2 weeks)

**Deliverables**:
1. Format specification (this document)
2. `rvids3 convert` command:
   - Input: video file + ELAN .eaf or Anvil .anvil
   - Output: RVIDS3 format
   - Video/audio segmentation
   - Index generation
3. `rvids3 stats` command
4. Unit tests

**Example**:
```bash
rvids3 convert session_video.mp4 session.eaf output.rvids3 \
  --segment-duration 30 \
  --video-codec h264 \
  --audio-codec opus
```

### Phase 2: Query Engine (2 weeks)

**Deliverables**:
1. `rvids3 query` command:
   - Time-based extraction
   - Annotation-based filtering
   - Speaker/participant filtering
   - Context inclusion (before/after)
2. `rvids3 extract` command:
   - Simple time-range extraction
   - Output to standard video formats
3. Integration tests

**Example**:
```bash
# Extract all reaching behaviors with 1-second context
rvids3 query session.rvids3 \
  --tier behavior \
  --code reach \
  --context 1.0 \
  --output clips/reaching/
```

### Phase 3: Export and Tool Integration (1 week)

**Deliverables**:
1. `rvids3 export` command:
   - ELAN .eaf export
   - Anvil .anvil export
   - CSV export (tabular)
   - JSON export (raw data)
2. Validation with real ELAN/Anvil projects
3. Documentation and examples

**Example**:
```bash
# Round-trip validation
rvids3 export session.rvids3 output.eaf --format elan
# Open in ELAN, verify all annotations present
```

### Phase 4: S3 Integration (1 week)

**Deliverables**:
1. S3 URI support for input/output
2. Selective segment downloads (range requests)
3. Streaming extraction
4. Cost analysis and benchmarks

**Example**:
```bash
# Convert and upload to S3
rvids3 convert video.mp4 annotations.eaf s3://bucket/session.rvids3

# Query directly from S3
rvids3 query s3://bucket/session.rvids3 \
  --tier behavior --code reach \
  --output /tmp/clips/
```

### Phase 5: Production Features (2 weeks)

**Deliverables**:
1. Batch processing for multiple sessions
2. Inter-rater reliability tools
3. Aggregate analysis across sessions
4. Python SDK for programmatic access
5. Web viewer (optional)

## Use Case Examples

### Example 1: Psychology Research - Parent-Child Interaction

**Dataset**: 50 sessions, 1 hour each, 100GB total

**Traditional Workflow**:
```bash
# Download all videos
aws s3 sync s3://study/videos/ ./videos/  # 100GB download
# Manually search in ELAN
# Extract clips manually
# Total time: Days
```

**RVIDS3 Workflow**:
```bash
# Query all joint attention episodes
rvids3 query s3://study/sessions/*.rvids3 \
  --tier interaction \
  --code joint_attention \
  --output analysis/joint_attention/ \
  --aggregate

# Result: 200 clips extracted, 2GB total downloaded
# Savings: 98GB, 98% reduction
# Time: Minutes
```

### Example 2: Linguistics - Conversation Analysis

**Dataset**: 20 dyadic conversations, 30 minutes each, 40GB total

**Research Question**: Analyze turn-taking patterns

**RVIDS3 Workflow**:
```bash
# Extract all speaker turns with timing
rvids3 export s3://corpus/*.rvids3 turns.csv \
  --format csv \
  --speakers all \
  --include-timing

# Analyze in R/Python without downloading videos
# Only download specific clips for examples

# Data transfer: 50MB (metadata only)
# Savings: 99.875%
```

### Example 3: Education Research - Classroom Observation

**Dataset**: 100 classroom sessions, 45 minutes each, 200GB total

**Research Question**: Teacher questioning strategies

**RVIDS3 Workflow**:
```bash
# Extract all teacher questions
rvids3 query s3://classrooms/*.rvids3 \
  --speaker teacher \
  --transcript-contains "?" \
  --context 5.0 \
  --output analysis/questions/

# Export annotations for inter-rater reliability
rvids3 export s3://classrooms/*.rvids3 icc_data.csv \
  --format csv \
  --tier questioning_strategy

# Data transfer: 5GB (2.5% of dataset)
# Analysis time: Hours instead of weeks
```

## Cost Analysis

### Storage Costs

**1-Hour Research Session**:
- Video (1080p, 30fps): 2GB
- Audio (stereo, 48kHz): 65MB
- Transcript (JSON, compressed): 200KB
- Annotations (3 tiers, compressed): 300KB
- Indices: 100KB
- **Total**: ~2.07GB

**Cost**: 2.07GB × $0.023/GB/month = $0.048/month/session

**Traditional** (separate files):
- Video: 2GB
- Transcript: 5MB (Word doc)
- ELAN file: 2MB (XML)
- Metadata: 1MB (Excel)
- **Total**: ~2.008GB

**Storage savings**: Minimal (~3% larger due to indices)
**Key benefit**: Unified format, not storage reduction

### Query Costs

**Scenario**: Extract 50 behavioral instances from 50 sessions (50 hours total)

**Traditional**:
- Download all videos: 100GB
- Cost: 100GB × $0.09/GB = $9.00
- Time: Hours

**RVIDS3**:
- Download only relevant segments: 2GB
- Cost: 2GB × $0.09/GB = $0.18
- Time: Minutes

**Savings**: $8.82 per query (98% reduction)

### Annual Research Project Costs

**Scenario**: 100 sessions, 1000 queries/year

**Traditional**:
- Storage: 200GB × $0.023/mo × 12 = $55.20
- Queries: 1000 × $9.00 = $9,000
- **Total**: $9,055.20

**RVIDS3**:
- Storage: 210GB × $0.023/mo × 12 = $57.96
- Queries: 1000 × $0.18 = $180
- **Total**: $237.96

**Annual savings**: $8,817.24 (97.4% reduction)

## Tool Ecosystem

### Core CLI Tool

**Commands**:
- `rvids3 convert` - Create RVIDS3 from video + annotations
- `rvids3 query` - Extract clips by annotations/time/speaker
- `rvids3 extract` - Simple time-based extraction
- `rvids3 export` - Convert to ELAN/Anvil/CSV
- `rvids3 stats` - Dataset statistics
- `rvids3 validate` - Verify format integrity

### Python SDK

```python
from rvids3 import Session

# Open session (local or S3)
session = Session("s3://bucket/session.rvids3")

# Query by annotation
clips = session.query(tier="behavior", code="reach", context=1.0)
for clip in clips:
    print(f"Found at {clip.start_time}s: {clip.metadata}")
    clip.save(f"clip_{clip.id}.mp4")

# Get transcript
transcript = session.transcript(speaker="child_01")
for utterance in transcript:
    print(f"[{utterance.start_time}s] {utterance.text}")

# Export annotations
session.export_tier("behavior", "behaviors.csv", format="csv")
```

### Integration with Existing Tools

**ELAN**: Import/export via XML conversion
**Anvil**: Import/export via XML conversion
**Praat**: Export TextGrid for prosody analysis
**CHILDES**: Export CHAT format for language analysis

## Technical Specifications

### Video Encoding

- **Container**: MP4 (H.264 + AAC/Opus)
- **Codec**: H.264 High Profile, Level 4.0
- **Resolution**: Preserve original (typical: 1080p, 720p, 480p)
- **Frame Rate**: Preserve original (typical: 30fps, 25fps)
- **Bitrate**: 2-4 Mbps (1080p), 1-2 Mbps (720p)
- **Keyframe Interval**: 1 second (30 frames at 30fps)

### Audio Encoding

- **Container**: Opus in Ogg or standalone
- **Codec**: Opus
- **Bitrate**: 32-64 kbps (mono), 64-96 kbps (stereo)
- **Sample Rate**: 48 kHz (Opus native)

### Compression

- **Annotations**: zstd level 3 (fast decode, good ratio)
- **Indices**: Uncompressed binary (fast access)
- **Metadata**: Uncompressed JSON (small files)

## Limitations and Considerations

### Format Limitations

1. **Segment Granularity**: Minimum extraction unit is ~1 second (keyframe interval)
2. **Video Quality**: Re-encoding may slightly reduce quality (use high bitrate to minimize)
3. **File Size**: Indices add ~5MB overhead per hour of video

### Research Workflow Considerations

1. **Annotation in ELAN/Anvil First**: Convert to RVIDS3 after annotation
2. **Inter-Rater Reliability**: Use RVIDS3 export to standard formats for ICC calculation
3. **Backup Original Files**: RVIDS3 is for analysis/sharing, keep originals

### Privacy and Ethics

1. **Identifiable Faces**: Video contains identifiable information
2. **Data Sharing**: Obtain proper IRB approval before sharing RVIDS3 files
3. **Cloud Storage**: Ensure compliance with institutional policies for sensitive data
4. **Access Control**: Use S3 bucket policies to restrict access

## Comparison with Alternatives

### RVIDS3 vs. ELAN

| Feature | RVIDS3 | ELAN |
|---------|--------|------|
| Format | Unified (video + annotations) | Separate (XML + video) |
| Cloud-native | Yes (S3 direct) | No |
| Selective access | Yes (annotation-based) | No (download all) |
| Programmatic API | Yes (Python, CLI) | Limited (XML parsing) |
| Tool support | Export to ELAN | Native |
| Sharing | Single file/directory | Multiple files |

**Use Case**: RVIDS3 for cloud-based analysis and sharing, ELAN for annotation

### RVIDS3 vs. Datavyu

| Feature | RVIDS3 | Datavyu |
|---------|--------|----------|
| Format | Cloud-native | Local files |
| Multi-tier | Yes | Yes |
| Reliability coding | Export to standard tools | Built-in |
| Scripting | Python SDK | Ruby API |
| Cost | S3 costs | Free software |

**Use Case**: RVIDS3 for large-scale studies, Datavyu for small local projects

## Future Enhancements

### Phase 6 (Future)

1. **Real-Time Annotation Sync**: Multi-coder real-time collaboration
2. **Machine Learning Integration**: Automatic behavior detection, transcript generation
3. **Cross-Session Queries**: "Find behavior X across all sessions in study"
4. **Web-Based Viewer**: Browser-based playback and annotation
5. **Multi-Camera Support**: Synchronized multi-angle video
6. **Physiological Data**: Heart rate, EDA, eye-tracking integration

## Summary

RVIDS3 provides research-focused video management with:

- **Unified format** for video, audio, transcripts, and annotations
- **Cloud-native architecture** with S3 direct access
- **Annotation-based queries** saving 98% on data transfer
- **Tool compatibility** with ELAN, Anvil, and standard formats
- **Cost savings** of $8,000+ per year for typical research projects
- **Reproducible workflows** with CLI and Python SDK

**Target Users**: Psychology researchers, linguists, education researchers, behavioral scientists

**Status**: Design complete, ready for implementation
