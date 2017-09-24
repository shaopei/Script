f = 'TssPair.d150_dRegion250_withStrandSpecific.GroseqDiscordanceInterestingHetsNoWeird.bed'

with open(f[:-1]+'modified.bed', 'w') as out:
    for l in open(f,'U').readlines():
        ll = l.strip().split('\t')
        out.write('\t'.join(ll[0:3]))
        out.write('\n')