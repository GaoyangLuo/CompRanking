#!/bin/bash


while getopts "p:i" option; do
	case "${option}" in
		p) PREFIX=${OPTARG};;
		i) INPUT_DIR=${OPTARG};;
		*) exit 1;;
	esac
done


workdir=${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file

# 遍历所有 .faa 文件
for i in "$workdir"/*.fna2faa.faa; do
    filename=$(basename "$i")         # 获取文件名
    base=${filename%%.fna2faa.faa}    # 去掉 .fna2faa.faa 后缀

    # 生成 ID 映射文件
    awk -F'[>_ ]' '
    /^>/ {
        contig = $2"_"$3;
        if (!(contig in count)) count[contig] = 1;
        else count[contig]++;
        old_id = $2"_"$3"_"$4"_"$5;
        new_id = contig"_"count[contig];
        print old_id, new_id;
    }' "$i" > "$workdir/id_map_${base}.txt"

    # 用 awk 进行整体替换，避免多次 sed 调用
    awk -v mapfile="$workdir/id_map_${base}.txt" '
    BEGIN {
        while ((getline < mapfile) > 0) id_map[$1] = $2;
        close(mapfile);
    }
    {
        if ($0 ~ /^>/) {
            split($0, arr, " ");
            old_id = substr(arr[1], 2);  # 去掉 ">"
            if (old_id in id_map) {
                arr[1] = ">" id_map[old_id];  # 替换 ID
            }
            $0 = arr[1];
            for (i=2; i<=length(arr); i++) $0 = $0 " " arr[i]; 
        }
        print;
    }' "$i" > "$workdir/${base}.faa"

done