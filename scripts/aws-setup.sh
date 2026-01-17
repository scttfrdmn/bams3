#!/bin/bash
set -e

# BAMS3 WGS Testing - AWS Infrastructure Setup
# Region: us-west-2
# Profile: aws

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║    BAMS3 WGS Testing - AWS Infrastructure Setup               ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Configuration
PROFILE="aws"
REGION="us-east-1"
BUCKET_NAME="bams3-testing-${USER}"
INSTANCE_TYPE="c5.9xlarge"
VOLUME_SIZE=500

echo -e "${BLUE}Configuration:${NC}"
echo "  AWS Profile:   $PROFILE"
echo "  Region:        $REGION (same as 1000 Genomes data)"
echo "  S3 Bucket:     $BUCKET_NAME"
echo "  Instance Type: $INSTANCE_TYPE"
echo "  Disk Size:     ${VOLUME_SIZE}GB"
echo ""

# Check for required parameters
if [ -z "$1" ] || [ -z "$2" ]; then
    echo -e "${RED}Error: Missing required parameters${NC}"
    echo ""
    echo "Usage: $0 <key-name> <security-group-id>"
    echo ""
    echo "Example:"
    echo "  $0 my-ec2-key sg-0123456789abcdef0"
    echo ""
    echo "Notes:"
    echo "  - Key pair must exist in us-west-2"
    echo "  - Security group must allow SSH (port 22)"
    echo "  - IAM role 'BAMS3-Testing-Role' must exist (see plan for permissions)"
    exit 1
fi

KEY_NAME="$1"
SG_ID="$2"

# Step 1: Create S3 bucket
echo -e "${BLUE}[1/3] Creating S3 bucket: $BUCKET_NAME${NC}"
if aws --profile $PROFILE s3 ls "s3://$BUCKET_NAME" 2>/dev/null; then
    echo -e "${YELLOW}  Bucket already exists, skipping${NC}"
else
    aws --profile $PROFILE s3 mb "s3://$BUCKET_NAME" --region $REGION
    echo -e "${GREEN}✓${NC} Bucket created"
fi
echo ""

# Step 2: Verify IAM role exists
echo -e "${BLUE}[2/3] Verifying IAM role: BAMS3-Testing-Role${NC}"
if aws --profile $PROFILE iam get-role --role-name BAMS3-Testing-Role >/dev/null 2>&1; then
    echo -e "${GREEN}✓${NC} IAM role exists"
else
    echo -e "${RED}✗ IAM role 'BAMS3-Testing-Role' not found${NC}"
    echo ""
    echo "Please create the role with the following permissions:"
    echo "  - S3 read from s3://1000genomes/*"
    echo "  - S3 read/write to s3://bams3-testing-*/*"
    echo "  - CloudWatch logs write"
    echo ""
    exit 1
fi
echo ""

# Step 3: Launch EC2 instance
echo -e "${BLUE}[3/3] Launching EC2 instance${NC}"
echo "  Instance type: $INSTANCE_TYPE"
echo "  Key pair:      $KEY_NAME"
echo "  Security group: $SG_ID"
echo ""

INSTANCE_ID=$(aws --profile $PROFILE ec2 run-instances \
  --image-id ami-0c55b159cbfafe1f0 \
  --instance-type $INSTANCE_TYPE \
  --region $REGION \
  --key-name "$KEY_NAME" \
  --security-group-ids "$SG_ID" \
  --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":$VOLUME_SIZE,\"VolumeType\":\"gp3\"}}]" \
  --iam-instance-profile Name=BAMS3-Testing-Role \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=BAMS3-WGS-Testing}]" \
  --query 'Instances[0].InstanceId' \
  --output text)

echo -e "${GREEN}✓${NC} Instance launched: $INSTANCE_ID"
echo ""

# Wait for instance to be running
echo -e "${BLUE}Waiting for instance to start...${NC}"
aws --profile $PROFILE ec2 wait instance-running --instance-ids "$INSTANCE_ID" --region $REGION
echo -e "${GREEN}✓${NC} Instance running"
echo ""

# Get public IP
PUBLIC_IP=$(aws --profile $PROFILE ec2 describe-instances \
  --instance-ids "$INSTANCE_ID" \
  --region $REGION \
  --query 'Reservations[0].Instances[0].PublicIpAddress' \
  --output text)

echo -e "${BOLD}${GREEN}Setup Complete!${NC}"
echo ""
echo -e "${BLUE}Instance Details:${NC}"
echo "  Instance ID: $INSTANCE_ID"
echo "  Public IP:   $PUBLIC_IP"
echo "  Region:      $REGION"
echo ""
echo -e "${BLUE}Next Steps:${NC}"
echo "  1. SSH to instance:"
echo "     ssh -i ~/.ssh/${KEY_NAME}.pem ec2-user@${PUBLIC_IP}"
echo ""
echo "  2. Copy and run ec2-setup-userdata.sh on the instance:"
echo "     scp -i ~/.ssh/${KEY_NAME}.pem scripts/ec2-setup-userdata.sh ec2-user@${PUBLIC_IP}:~/"
echo "     ssh -i ~/.ssh/${KEY_NAME}.pem ec2-user@${PUBLIC_IP}"
echo "     chmod +x ec2-setup-userdata.sh"
echo "     sudo ./ec2-setup-userdata.sh"
echo ""
echo "  3. Run tests (once setup is complete):"
echo "     ./run-chr22-test.sh"
echo "     ./run-full-wgs-test.sh"
echo ""
echo -e "${YELLOW}Estimated setup time: ~90 minutes (reference indexing)${NC}"
echo -e "${YELLOW}Estimated cost: ~$11 for full test suite${NC}"
echo ""

# Save instance info
cat > instance-info.txt << EOF
Instance ID: $INSTANCE_ID
Public IP: $PUBLIC_IP
Region: $REGION
Bucket: $BUCKET_NAME
Launched: $(date)
EOF

echo -e "${GREEN}✓${NC} Instance info saved to instance-info.txt"
echo ""
