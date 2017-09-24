f = 'maternal.chain.header.txt'

for l in open(f).readlines():
    ll = l.strip().split(' ')
    print 'sed -i \'s/'+' '.join(ll[0:3])+'/'+' '.join(ll[0:2])+' chr'+ll[2]+'/g\' maternal.chain'
    
    #sed -i 's/chain 102506170 15/chain 102506170 chr15/g' paternal.chain