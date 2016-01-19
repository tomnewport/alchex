def martini_type_similarity(type1, type2):
    if type1 == type2:
        return 1
    else:
        differences = {
            "C1.C3" : 0.6,
            "P4.P5" : 0.6
        }
        return getattr(differences, ".".join(sorted([type1, type2])), 0)