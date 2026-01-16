# GitHub Issues for BAMS3 Project

This directory contains prepared GitHub issue templates for remaining work.

## Creating Issues on GitHub

Once the repository is pushed to GitHub, create issues using these templates:

### Option 1: Manual Creation
1. Go to your repository on GitHub
2. Click "Issues" → "New Issue"
3. Copy content from each `.md` file
4. Add appropriate labels and assignees

### Option 2: Using GitHub CLI
```bash
# Install GitHub CLI if needed: brew install gh

# Authenticate
gh auth login

# Create issues from templates
gh issue create --title "Complete disk spill implementation" --body-file .github/issues/01-disk-spill.md --label enhancement,performance
gh issue create --title "End-to-end testing with BWA aligner" --body-file .github/issues/02-bwa-testing.md --label testing,priority-high
gh issue create --title "S3 direct integration" --body-file .github/issues/03-s3-integration.md --label enhancement,cloud,priority-high
gh issue create --title "Create Nextflow reference pipeline" --body-file .github/issues/04-nextflow-pipeline.md --label pipeline,documentation,priority-high
gh issue create --title "Implement BAM export" --body-file .github/issues/05-bam-export.md --label enhancement,integration,priority-high
gh issue create --title "Performance benchmarking" --body-file .github/issues/06-benchmarking.md --label performance,documentation,priority-high
```

## Issue Summary

| # | Title | Priority | Status |
|---|-------|----------|--------|
| 1 | Complete disk spill implementation | Medium | Infrastructure ready |
| 2 | End-to-end testing with BWA aligner | High | Blocked on BWA install |
| 3 | S3 direct integration | High | Not started |
| 4 | Create Nextflow reference pipeline | High | Not started |
| 5 | Implement BAM export | High | Not started |
| 6 | Performance benchmarking | High | Not started |

## Suggested Labels

Create these labels in your repository:
- `enhancement` - New feature or improvement
- `testing` - Testing-related tasks
- `performance` - Performance optimization
- `cloud` - Cloud integration
- `pipeline` - Workflow/pipeline related
- `documentation` - Documentation improvements
- `integration` - Tool integration
- `priority-high` - High priority
- `priority-medium` - Medium priority
- `priority-low` - Low priority

## Project Board Setup

Consider creating a GitHub Project board with columns:
- **Backlog** - Not yet started
- **In Progress** - Currently being worked on
- **Blocked** - Waiting on dependencies
- **Review** - Ready for review
- **Done** - Completed

## Notes

These issues represent the core remaining work to complete the BAMS3 demonstration project. They are organized by priority and include detailed acceptance criteria.
