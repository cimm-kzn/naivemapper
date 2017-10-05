# from CGRtools.strings import hash_cgr_string, get_morgan, get_cgr_string
from CGRtools.containers import MoleculeContainer
from CGRtools.files.RDFrw import RDFread, RDFwrite
from CGRtools.preparer import CGRcombo
import sys

cgr = CGRcombo()

with open(sys.argv[1], encoding='cp1251') as fr, open(sys.argv[2], encoding='cp1251') as fw:
    ok, nok = 0, 0
    er = []

    for i, (test, pred) in enumerate(zip(RDFread(fr), RDFread(fw)), start=1):
        predHash = cgr.getCGR(pred).get_fear_hash()
        testHash = cgr.getCGR(test).get_fear_hash()
        if predHash == testHash:
            ok += 1
        else:
            nok += 1
            er.append(i)

    print("Percentage\n\tO'k: %0.5f , \nNot O'k: %0.5f" % ((ok*100/(ok + nok)), (nok*100/(ok + nok))))
    print(len(er), '\n', er)
"""
with open(sys.argv[1], encoding='cp1251') as fr, open(sys.argv[2], encoding='cp1251') as fw:
    l = [5, 10, 12, 14, 15, 18, 21, 22, 27, 30, 31, 32, 33, 37, 38, 39, 43, 45, 48, 49, 50, 51, 56, 57, 65, 68, 75, 78, 79, 81, 82, 84, 85, 92, 94, 99, 101, 102, 106, 107, 115, 117, 118, 119, 128, 130, 131, 135, 136, 145, 148, 152, 153, 155, 158, 159, 161, 166, 170, 173, 178, 179, 185, 189, 191, 192, 197, 199, 202, 208, 209, 210, 211, 215, 216, 219, 220, 221, 222, 223, 226, 229, 231, 233, 236, 238, 240, 242, 243, 246, 248, 249, 251, 255, 260, 261, 262, 265, 269, 270, 271, 275, 276, 277, 282, 283, 285, 287, 289, 296, 297, 301, 305, 307, 309, 310, 317, 318, 324, 325, 326, 330, 333, 337, 341, 342, 343, 345, 346, 348, 349, 351, 354, 355, 358, 362, 364, 366, 369, 377, 379, 380, 381, 384, 392, 394, 396, 398, 399, 406, 408, 409, 417, 421, 427, 431, 433, 435, 436, 438, 440, 442, 445, 448, 449, 450, 451, 452, 453, 455, 456, 462, 468, 469, 470, 476, 479, 480, 484, 486, 488, 489, 490, 491, 492, 494, 495, 499, 502, 504, 505, 506, 509, 514, 516, 517, 519, 524, 527, 532, 534, 537, 540, 544, 545, 547, 549, 554, 555, 558, 566, 568, 569, 575, 577, 579, 583, 585, 586, 592, 597, 599, 603, 615, 618, 619, 624, 629, 630, 632, 635, 642, 647, 652, 655, 656, 673, 674, 678, 680, 682, 683, 686, 689, 692, 693, 702, 703, 704, 707, 708, 709, 710, 713, 714, 716, 718, 720, 721, 725, 726, 727, 729, 733, 737, 738, 742, 743, 744, 747, 749, 754, 760, 762, 763, 765, 767, 771, 776, 777, 780, 781, 782, 791, 793, 797, 798, 801, 803, 805, 807, 808, 810, 811, 814, 816, 820, 822, 823, 825, 827, 828, 829, 830, 832, 836, 837, 839, 840, 845, 847, 850, 857, 858, 859, 860, 861, 863, 871, 874, 881, 882, 883, 887, 888, 891, 895, 903, 904, 905, 907, 908, 909, 915]
    ii = 0
    er2 = []
    predFile = RDFread(fw)
    ok, nok = 0, 0
    er = []
    for i, test in enumerate(RDFread(fr), start=1):  # (test, pred) in enumerate(zip(RDFread(fr), RDFread(fw)), start=1):
        if i not in l:
            pred = next(predFile)
            ii += 1
            predHash = cgr.getCGR(pred).get_fear_hash()
            testHash = cgr.getCGR(test).get_fear_hash()
            if predHash == testHash:
                ok += 1
            else:
                nok += 1
                er.append(i)
                er2.append(ii)
        else:
            nok += 1
            er.append(i)

    print("Percentage\n\tO'k: %0.5f , \nNot O'k: %0.5f" % ((ok*100/(ok + nok)), (nok*100/(ok + nok))))
    print(len(er), '\n', er)
    print(ii, '\t', len(er2), '\n', er2)
"""
"""
with open("test/Etherification/esterification.rdf", encoding='cp1251') as fr, open(input('Name file:'), "w") as fw:
    uniqueHash = {}
    print('Seeking unique items')
    for num, reaction in enumerate(RDFread(fr), start=1):
        rHash = cgr.getCGR(reaction).get_fear_hash()
        uniqueHash[rHash] = reaction

    print('Record file')
    outputdata = RDFwrite(fw)
    for v in uniqueHash.values():
        outputdata.write(v)

    print(len(uniqueHash), ' unique reactions')
"""
"""
with open('test/molecule_etr.rdf') as f, open('test/reac_etr.rdf') as fr, open('test/reac_etr_out.rdf', 'w') as fw:
    r = CGR.getCGR(RDFread(f).read()[0])
    rHash = fear.get_cgr_string(r.subgraph(fear.get_center_atoms(r)))
    out = RDFwrite(fw)
    for num, reac in enumerate(RDFread(fr).read(), start=1):
        reacCGR = CGR.getCGR(reac)
        reacCenter = reacCGR.subgraph(fear.get_center_atoms(reacCGR))
        try:
            rcHash = fear.get_cgr_string(reacCenter)
            if rHash in rcHash:
                out.write(reac)
                print(num)
                # print(rcHash)
        except Exception as e:
            print(e)
"""
