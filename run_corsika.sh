#!/bin/bash

# 1. รับค่าพารามิเตอร์
PART_NAME=$1
E_VAL=$2
ZENITH_VAL=$3

# 2. ทำการ Compile โปรแกรมบนเครื่อง Worker
# การรันที่นี่จะปลอดภัยเพราะแต่ละ Job มีโฟลเดอร์ชั่วคราว (Scratch Space) ของตัวเอง
echo "Compiling for $PART_NAME, $E_VAL, $ZENITH_VAL..."
make clean
make PARTICLE=$PART_NAME E_VAL=$E_VAL ZENITH=sin2_$ZENITH_VAL

# ตรวจสอบว่า compile สำเร็จหรือไม่
if [ ! -f "./corsikaReader" ]; then
    echo "Error: Compilation failed!"
    exit 1
fi

# 3. กำหนด Path ของข้อมูล (ซึ่งปกติจะอยู่ใน Network Storage ที่ทุกเครื่องมองเห็น)
BASE_PATH="/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/${PART_NAME}/${E_VAL}/sin2_${ZENITH_VAL}"

# 4. วนลูปการรันโปรแกรม
for file_path in ${BASE_PATH}/000???/DAT000???; do
    if [[ -f "$file_path" ]]; then
        echo "Running: ./corsikaReader $file_path --thinned"
        ./corsikaReader "$file_path" --thinned
    else
        echo "Warning: No files found in $BASE_PATH"
    fi
done

