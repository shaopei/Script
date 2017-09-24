from sys import argv

f = argv[1]


with open(f, 'U') as f_input:
    with open (f+'_plus', 'w') as p_out:
        with open (f+'_minus', 'w') as m_out:
            with open (f+'_question', 'w') as q_out:
                for l in f_input.readlines():
                    s = l.split('\t')[1]
                    print s 
                    if s == '+':
                        p_out.write(l)
                    elif s == '-':
                        m_out.write(l)
                    else:
                        q_out.write(l)
