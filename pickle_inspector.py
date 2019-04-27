def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

def clever_print(x):
    a = namestr(x, globals())[0]
    print(a, x)

import pickle
results = pickle.load(open("storage.p", 'rb'))



show_all = True
if show_all:
    c = 0
    for r in results:
        if r["blast"].q2hmap == r["sail"].q2hmap:
            print("---- " + str(c) + " ----- PERFECT q2hmap MATCH")
        else:
            print("---- " + str(c) + " -----  SEE q2hmap MIS-MATCH")
        c = c + 1

    c = 0
    for r in results:
        if not r["blast"].q2hmap == r["sail"].q2hmap:
            print("---- "+ str(c) + " -----")

            for i in ["blast", "sail"]:
                print("h_strand", r[i].h_strand)
                print("q_strand", r[i].q_strand)
                print("h_start", r[i].h_start)

                print(">>> " + i + " q2hmap >>> " ,
                      r[i].q2hmap.keys()[0],
                      r[i].q2hmap[r[i].q2hmap.keys()[0]],
                      r[i].q2hmap.keys()[-1],
                      r[i].q2hmap[r[i].q2hmap.keys()[-1]])
                print(r[i].q_align)
                print(r[i].middleseq)
                print(r[i].h_align)

            # print(">>> SAIL q2hmap >>> ",
            #       r["sail"].q2hmap.keys()[0],
            #       r["sail"].q2hmap[r["sail"].q2hmap.keys()[0]],
            #       r["sail"].q2hmap.keys()[-1],
            #       r["sail"].q2hmap[r["sail"].q2hmap.keys()[-1]])
            # print(r["sail"].q_align)
            # print(r["sail"].middleseq)
            # print(r["sail"].h_align)
            print
            print

        c = c + 1