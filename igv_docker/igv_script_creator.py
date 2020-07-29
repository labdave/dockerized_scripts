import argparse
import os

""" GET ARGUMENTS FROM COMMAND LINE"""

parser = argparse.ArgumentParser(description='FILE FORMAT REQUIRED: SV_CALLING_ID CHR1 POS1 CHR2 POS2')
parser.add_argument('-b', '--bam', required=True, help='bam file to get snapshots from')
parser.add_argument('-q', '--squished', action='store_true', help='allow squished output')
parser.add_argument('-s', '--split', action='store_true', help='allow split screen with mate')
parser.add_argument('-g', '--group', action='store_true', help='allow group by mate chromosome')
parser.add_argument('-f', '--file', required=True, help='translocation input file')
parser.add_argument('-z', '--zoom', required=True, help='zoom amount', type=int)
parser.add_argument('-d', '--dir', required=True, help='image locations')
parser.add_argument('-o', '--out', required=True, help='output name')
parser.add_argument('-r', '--repeat', required=True, help='repeat masker exclude list')
parser.add_argument('-t', '--target', required=True, help='target bed')
args = parser.parse_args()


""" PARSE INPUT FILE """

analyis_dict = dict()
with open(args.file, 'r') as f:
    for line in f:
        # expected line format:
        # SV_CALLING_ID <tab> CHR1 <tab> POS1 <tab> CHR2 <tab> POS2
        # 
        # example file
        # 32489 chr14   105864253   chr18   63126261
        # 32489 chr8    127795134   chr22   22893781
        # 32489 chr8    127795134   chr22   22895373
        # 32489 chr3    111555239   chr8    127521650
        # 32489 chr8    127919318   chr22   22881132
        # 32489 chr8    127919351   chr22   22881114
        # 32489 chr8    127919318   chr22   22881132
        # 32489 chr8    127919320   chr22   22881132
        # 32489 chr14   106113360   chr16   33841302
        # 32490 chr8    127736611   chr14   105859660
        # 32490 chr8    127736628   chr14   105859664
        # 32490 chr8    127736611   chr14   105859660
        # 32490 chr8    127736628   chr14   105859664
        # 32490 chr8    127736611   chr14   105859660
        # 32490 chr8    127736628   chr14   105859664
        # 32490 chr8    127736611   chr14   105859660
        # 32490 chr8    127736628   chr14   105859664
        # 32490 chr3    111555239   chr8    127521655
        # 32490 chr8    127733934   chr14   105864213
        # 32490 chr6    115311362   chr14   106782735
        if 'id' in line.lower() or 'pos' in line.lower():
            continue
        line_arr = line.strip().split()
        analyis_id = line_arr[0]
        chr1 = line_arr[1]
        pos1 = line_arr[2]
        pos1a = int(pos1) - args.zoom
        pos1b = int(pos1) + args.zoom
        chr2 = line_arr[3]
        pos2 = line_arr[4]
        pos2a = int(pos2) - args.zoom
        pos2b = int(pos2) + args.zoom
        
        if analyis_id in analyis_dict:
            analyis_dict[analyis_id].append([chr1, pos1, pos1a, pos1b, chr2, pos2, pos2a, pos2b])
        else:
            analyis_dict[analyis_id] = [[chr1, pos1, pos1a, pos1b, chr2, pos2, pos2a, pos2b]]


""" START CREATING SCRIPT """

script = ''

# iterate through records
for analyis_id in analyis_dict:
    # start igv
    script += 'new\n'
    
    # load genome
    script += 'genome hg38\n'
    
    # load file
    script += 'load {0}\n'.format(args.bam)
    
    # load bed file
    script += 'load {}\n'.format(args.target)
    script += 'load {}\n'.format(args.repeat)

    # set screenshot directory
    script += 'snapshotDirectory '
    script += args.dir
    script += '\n'
    
    for record in analyis_dict[analyis_id]:
        chr1 = record[0]
        pos1 = record[1]
        pos1a = record[2]
        pos1b = record[3]
        chr2 = record[4]
        pos2 = record[5]
        pos2a = record[6]
        pos2b = record[7]
        
        # add gotos
        if args.split:
            script += 'goto {0}:{1}-{2} {3}:{4}-{5}\n'.format(chr1, pos1a, pos1b, chr2, pos2a, pos2b)
        else:
            script += 'goto {0}:{1}-{2}\n'.format(chr1, pos1a, pos1b)
        
        # increase panel height
        script += 'maxPanelHeight 4000\n'

        # add sort position
        script += 'sort position\n'

        # add squished or collapsed
        if args.squished:
            script += 'squish\n'
            script += 'maxPanelHeight 4000\n'
        else:
            script += 'expand\n'
        
        # group my mate chromosome
        if args.group:
            script += 'group MATE_CHROMOSOME\n'
        
        # fix stuff
        script += 'collapse Gene\n'
        script += 'expand Twist_8MB_panel_with_ERCCs.maskPAR.bed\n'
        script += 'expand hg38_repeat_masker.sorted.bed\n'

        # get snapshot
        script += 'snapshot {0}_{1}_{2}-{3}_{4}'.format(analyis_id, chr1, pos1, chr2, pos2)
        if args.split:
            script += '_split'
        if args.squished:
            script += '_squished'
        if args.group:
            script += '_group'
        script += '.png\n'
    script += 'exit\n'


""" WRITE TO SCRIPT FILE"""

with open(args.out, 'w') as f:
    f.write(script)