#!/bin/bash
set -e

# Run all BAMS3/VCFS3 demos
# Generates consolidated report

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║          BAMS3/VCFS3 Complete Demo Suite                      ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check mode
INTERACTIVE=true
if [ "$1" = "--ci-mode" ] || [ "$1" = "--non-interactive" ]; then
    INTERACTIVE=false
fi

# Demo list
DEMOS=(
    "01-basic-workflow:Basic Workflow:2"
    "02-sharding-performance:S3 Sharding Performance:3"
    "03-vcfs3-selective-access:VCFS3 Selective Access:3"
    "04-cloud-vs-traditional:Cloud vs Traditional:3"
)

TOTAL_DEMOS=${#DEMOS[@]}
COMPLETED=0
FAILED=0

echo "Running $TOTAL_DEMOS demos..."
echo ""

mkdir -p results
REPORT_FILE="results/demo-report-$(date +%Y%m%d-%H%M%S).txt"

{
    echo "BAMS3/VCFS3 Demo Suite Report"
    echo "Generated: $(date)"
    echo "======================================"
    echo ""
} > "$REPORT_FILE"

for demo_info in "${DEMOS[@]}"; do
    IFS=':' read -r demo_dir demo_name demo_time <<< "$demo_info"

    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}Demo: $demo_name${NC}"
    echo -e "Expected time: ~${demo_time} minutes"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""

    if [ "$INTERACTIVE" = true ]; then
        read -p "Run this demo? (y/n/q) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Qq]$ ]]; then
            echo "Quitting demo suite..."
            break
        elif [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo -e "${YELLOW}Skipped${NC}"
            echo ""
            continue
        fi
    fi

    # Run demo
    START=$(date +%s)
    if cd "$demo_dir" && ./run-demo.sh; then
        END=$(date +%s)
        ELAPSED=$((END - START))

        echo -e "${GREEN}✓ $demo_name completed in ${ELAPSED}s${NC}"
        COMPLETED=$((COMPLETED + 1))

        {
            echo "✓ $demo_name"
            echo "  Time: ${ELAPSED}s"
            echo ""
        } >> "../$REPORT_FILE"

        cd ..
    else
        echo -e "${RED}✗ $demo_name failed${NC}"
        FAILED=$((FAILED + 1))

        {
            echo "✗ $demo_name FAILED"
            echo ""
        } >> "../$REPORT_FILE"

        cd ..
    fi

    echo ""
    echo ""
done

# Summary
echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║                     Demo Suite Summary                        ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Total demos:    $TOTAL_DEMOS"
echo -e "Completed:      ${GREEN}$COMPLETED${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "Failed:         ${RED}$FAILED${NC}"
fi
echo ""
echo "Report saved to: $REPORT_FILE"
echo ""

# Generate HTML report (optional)
if command -v pandoc >/dev/null 2>&1; then
    HTML_REPORT="${REPORT_FILE%.txt}.html"
    pandoc "$REPORT_FILE" -o "$HTML_REPORT" 2>/dev/null || true
    if [ -f "$HTML_REPORT" ]; then
        echo "HTML report: $HTML_REPORT"
    fi
fi

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}${BOLD}All demos completed successfully!${NC}"
    exit 0
else
    echo -e "${RED}${BOLD}Some demos failed. Check the report for details.${NC}"
    exit 1
fi
