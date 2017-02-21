from itertools import product

f_txt = input('File name:')
with open(f_txt) as f:
    d1 = {'correct': [], 'incorrect': []}  # для корр./некорр. отображений запоминаем значение meanLikelihood
    d2 = {}  # для каждого meanLikelihood запоминаем № реакции и значение маппинга(корр./некорр)
    for line in f:
        l = line.split()
        i, m, r = int(l[0]), l[1], float(l[2])  # № реакции, значение маппинга(корр./некорр), значение meanLikelihood
        if r not in d2:
            d2[r] = [[], []]
        d1[m].append(r)  # для данного значения маппинга(корр./некорр), запоминаем значение meanLikelihood
        d2[r][0].append(i)  # для данного meanLikelihood запоминаем № реакции
        d2[r][1].append(m)  # для данного meanLikelihood запоминаем значение маппинга(корр./некорр)

    l = [(i, j) for i, j in product(d1['correct'], d1['incorrect'])]
    # все паросочетания meanLikelihood между корр. и некорр. значенями маппинга
    iap = len([1 for i, j in l if i < j])/len(l)  # IAP - Invariant Accuracy of Prediction

    tpr, fpr, bac, s1, s2 = 0, 0, {}, [], []
    for i in sorted(d2):  # в порядке возрастания значения meanLikelihood
        l1, l2 = d2[i][0], d2[i][1]  # списки реакций и значений маппинга(корр./некорр) для данного meanLikelihood

        s1.extend([l1[j] for j, k in enumerate(l2) if k == 'correct'])
        # добавление № реакций с корр.маппингом для данного meanLikelihood
        s2.extend([l1[j] for j, k in enumerate(l2) if k == 'incorrect'])
        # добавление № реакций с некорр.маппингом для данного meanLikelihood

        tpr += l2.count('correct')/len(d1['correct'])
        # доля корр.маппингов при данном meanLikelihood в общем количестве корр.маппов
        fpr += l2.count('incorrect')/len(d1['incorrect'])
        # доля некорр.маппингов при данном meanLikelihood в общем количестве некорр.маппов

        b = round(0.5*(tpr+(1-fpr)), 3)  # BAC - balanced accuracy
        if b not in bac:
            bac[b] = [i, sorted(s1.copy()), sorted(s2.copy())]
            # для данного BAC запоминаем meanLikelihood и реакции с корр./некорр. маппом (имеющих <=meanLikelihood)

    opt_b = max(bac)  # максимальное значение BAC
    tp, fp = bac[opt_b][1], bac[opt_b][2]  # списки р-ий с корр./некорр. маппом при opt(BAC)
    all_reac = set(s1).union(set(s2))  # список всех реакций
    fn, tn = set(s1).difference(set(tp)), set(s2).difference(set(fp))  # списки р-ий с корр./некорр. маппом вне opt(BAC)

    print('Model has: \nОшибок маппирования\t%d\nIAP\t%0.3f\nOptimal BAC\t%0.3f\nOptimal Likelihood\t%0.9f' %
          (len(s2), iap, opt_b, bac[opt_b][0]))
    print('TP\t%0.3f\n\t(%d / %d)' % (100*len(tp)/len(all_reac), len(tp), len(all_reac)))
    print('FP\t%0.3f\n\t(%d / %d)' % (100*len(fp)/len(all_reac), len(fp), len(all_reac)))
    print('TN\t%0.3f\n\t(%d / %d)' % (100*len(tn)/len(all_reac), len(tn), len(all_reac)))
    print('FN\t%0.3f\n\t(%d / %d)' % (100*len(fn)/len(all_reac), len(fn), len(all_reac)))
