#!/bin/bash


set -uo pipefail


WORKDIR="WORKDIR"

MD5_FILE=$WORKDIR/md5_igsr.txt
MD5_CHECK_RESULT="MD5_CHECK_RESULT.txt"
URL_LIST=$WORKDIR/file_list.txt

CPU=12
THREADS=8


cd $WORKDIR

mkdir -p logs
touch $MD5_CHECK_RESULT

########################################
# Step1 

extract_ok_files() {

    echo
    echo "Extracting OK files..."

    grep ": OK" $MD5_CHECK_RESULT | \
    awk -F: '{print $1}' > ok_files.txt

    OK_NUM=$(wc -l ok_files.txt | awk '{print $1}')

    echo "OK file number: $OK_NUM"
}

########################################
# Step2 

generate_failed_list() {

    echo
    echo "Generating failed download list..."

    > failed_download_list.txt

    while read -r url
    do
        file=$(basename "$url")

        if ! grep -Fxq "$file" ok_files.txt; then
            echo "$url" >> failed_download_list.txt
        fi

    done < $URL_LIST

    NEED=$(wc -l failed_download_list.txt | awk '{print $1}')

    echo "Need download: $NEED"
}

########################################
# Step3

download_one() {

    url=$1
    file=$(basename "$url")

    echo "Downloading $file"

    aria2c \
        -x $THREADS \
        -s $THREADS \
        -k 4M \
        -c \
        --file-allocation=none \
        --console-log-level=warn \
        --retry-wait=10 \
        --max-tries=0 \
        "$url" \
        > logs/${file}.log 2>&1

    md5_line=$(grep "$file" $MD5_FILE)

    if [ -z "$md5_line" ]; then

        echo "No md5 info: $file"
        echo "$url" >> download_failed_again.txt
        return

    fi

    echo "$md5_line" | md5sum -c - > logs/${file}.md5 2>&1

    if grep -q "OK" logs/${file}.md5; then

        echo "Success: $file"

        echo "$file: OK" >> $MD5_CHECK_RESULT

    else

        echo "MD5 failed: $file"

        rm -f "$file"

        echo "$url" >> download_failed_again.txt

    fi
}

export -f download_one
export THREADS
export MD5_FILE
export MD5_CHECK_RESULT

########################################
# Main

ROUND=1

while true
do

    echo
    echo "======================================"
    echo "ROUND $ROUND"
    echo "======================================"

    extract_ok_files

    generate_failed_list

    NEED=$(wc -l failed_download_list.txt | awk '{print $1}')

    if [ "$NEED" -eq 0 ]; then

        echo
        echo "######################################"
        echo "ALL FILES DOWNLOADED AND MD5 OK"
        echo "######################################"

        break

    fi

    > download_failed_again.txt

    echo
    echo "Parallel downloading..."

    cat failed_download_list.txt | \
    xargs -n 1 -P $CPU -I {} bash -c 'download_one "$@"' _ {}

    if [ -s download_failed_again.txt ]; then

        echo
        echo "Retry failed files in next round..."

    fi

    ROUND=$((ROUND+1))

done

echo
echo "======================================"
echo "DOWNLOAD COMPLETED"
echo "======================================"