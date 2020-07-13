# python igv_creator_test.py -f cc_temp -z 300 -d /Users/Sara/Documents/DDB/Validation/translocation/igv_script_creator/test -o test.script -q -s -g
import argparse
import os

""" GET ARGUMENTS FROM COMMAND LINE"""

parser = argparse.ArgumentParser(description='FILE FORMAT REQUIRED: SV_CALLING_ID CHR1 POS1 CHR2 POS2 FISH_CALL')
parser.add_argument('-q', '--squished', action='store_true', help='allow squished output')
parser.add_argument('-s', '--split', action='store_true', help='allow split screen with mate')
parser.add_argument('-g', '--group', action='store_true', help='allow group by mate chromosome')
parser.add_argument('-f', '--file', required=True, help='translocation input file')
parser.add_argument('-z', '--zoom', required=True, help='zoom amount', type=int)
parser.add_argument('-d', '--dir', required=True, help='image locations')
parser.add_argument('-o', '--out', required=True, help='output name')
args = parser.parse_args()


""" PARSE INPUT FILE """

analysis_dict = dict()
with open(args.file, 'r') as f:
    for line in f:
        # expected line format:
        # SV_CALLING_ID <tab> CHR1 <tab> POS1 <tab> CHR2 <tab> POS2 <tab> FISH_CALL
        #
        # example file
        # 32489 chr14   105864253   chr18   63126261	BCL2
        # 32489 chr8    127795134   chr22   22893781	MYC
        # 32489 chr8    127919318   chr22   22881132	MYC
        # 32489 chr8    127919320   chr22   22881132	MYC
        # 32490 chr8    127736611   chr14   105859660	MYC
        # 32490 chr8    127736628   chr14   105859664	MYC
        # 32490 chr8    127736611   chr14   105859660	MYC
        # 32490 chr8    127736628   chr14   105859664	MYC
       
        if 'id' in line.lower() or 'pos' in line.lower():
            continue
        line_arr = line.strip().split()
        analysis_id = line_arr[0]
        chr1 = line_arr[1]
        pos1 = line_arr[2]
        pos1a = int(pos1) - args.zoom
        pos1b = int(pos1) + args.zoom
        chr2 = line_arr[3]
        pos2 = line_arr[4]
        pos2a = int(pos2) - args.zoom
        pos2b = int(pos2) + args.zoom
        fish = line_arr[5]

        if analysis_id in analysis_dict:
            analysis_dict[analysis_id].append([chr1, pos1, pos1a, pos1b, chr2, pos2, pos2a, pos2b,fish])
        else:
            analysis_dict[analysis_id] = [[chr1, pos1, pos1a, pos1b, chr2, pos2, pos2a, pos2b, fish]]


""" START CREATING SCRIPT """

script = ''

# home path
home ='s3://davelab-analysis-results/'

# iterate through records
for analysis_id in analysis_dict:
    # start igv
    script += 'new\n'

    # load genome
    script += 'genome hg38\n'

    print(analysis_id)

    # grep to get exact file path because aws does not allow wildcards
    cmd = 'aws s3 ls {0}{1}/ --recursive | grep discordant_reads | grep -v bai | cut -c 32-'.format(home, analysis_id)

    # get exact file name
    file_path = home + os.popen(cmd).read().strip()

    # load file
    script += 'load {0} index={1}\n'.format(file_path, file_path+'.bai')

    # set screenshot directory
    script += 'snapshotDirectory '
    script += args.dir
    script += '/{0}'.format(analysis_id)
    script += '\n'
    
    # make folders
    os.system('mkdir -p {0}/{1}'.format(args.dir, analysis_id))

    # load bed file
    script += 'load gs://davelab_data/ref/human/hg38/capture_baits/twist/maskPAR/Twist_8MB_panel_with_ERCCs.maskPAR.bed\n'
    script += 'load gs://davelab_data/ref/naresh_sv/hg38_repeat_masker.sorted.bed\n'
    script += 'load gs://davelab_data/ref/naresh_sv/hg38_segmental_dups.sorted.bed\n'

    for record in analysis_dict[analysis_id]:
        chr1 = record[0]
        pos1 = record[1]
        pos1a = record[2]
        pos1b = record[3]
        chr2 = record[4]
        pos2 = record[5]
        pos2a = record[6]
        pos2b = record[7]
        fish = record[8]

        # add gotos
        if args.split:
            script += 'goto {0}:{1}-{2} {3}:{4}-{5}\n'.format(chr1, pos1a, pos1b, chr2, pos2a, pos2b)
        else:
            script += 'goto {0}:{1}-{2}\n'.format(chr1, pos1a, pos1b)

        # add sort position
        script += 'sort position\n'

        # add squished or collapsed
        if args.squished:
            script += 'squish\n'
        else:
            script += 'collapse\n'
            
        #expand annotations
        script += 'expand ref%2Fnaresh_sv%2Fhg38_repeat_masker.sorted.bed\n'
        script += 'squish ref%2Fnaresh_sv%2Fhg38_segmental_dups.sorted.bed\n'
        script += 'expand ref%2Fhuman%2Fhg38%2Fcapture_baits%2Ftwist%2FmaskPAR%2FTwist_8MB_panel_with_ERCCs.maskPAR.bed\n'
        script += 'collapse Gene\n'

        # group my mate chromosome
        if args.group:
            script += 'group MATE_CHROMOSOME\n'

        # get snapshot
        script += 'snapshot {0}_{1}_{2}_{3}-{4}_{5}'.format(fish, analysis_id, chr1, pos1, chr2, pos2)
        if args.split:
            script += '_split'
        if args.squished:
            script += '_squished'
        if args.group:
            script += '_group'
        script += '.png\n'



""" WRITE TO SCRIPT FILE"""

with open(args.out, 'w') as f:
    f.write(script)
