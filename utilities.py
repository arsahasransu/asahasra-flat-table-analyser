import time


# Decorator to measure the execution time of a function
def time_eval(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Execution time: {execution_time:.6f} seconds")
        return result
    return wrapper


threep_vars = {
    'pt': (2001, -1, 2000),
    'eta': (54, -2.7, 2.7),
    'phi': (66, -3.3, 3.3),
}


charged_threep_vars = threep_vars | {
    'charge': (5, -2, 3),
    'vz': (300, -15, 15)
}


genel_vars = charged_threep_vars | {
    'caloeta': (54, -2.7, 2.7),
    'calophi': (66, -3.3, 3.3),
    'prompt': (5, -1, 4)
}


tkell2_vars = charged_threep_vars | {
    'hwQual': (10, -1, 9),
    'tkPt': (101, -1, 100),
    'caloEta': (54, -2.7, 2.7),
    'tkEta': (54, -2.7, 2.7),
    'caloPhi': (66, -3.3, 3.3),
    'tkPhi': (66, -3.3, 3.3),
    'pfIso': (2000, 0, 2),
    'puppiIso': (10000, 0, 10),
    'tkIso': (10000, 0, 10)
}


puppi_vars = threep_vars | {
    'mass': (500, 0, 0.5),
    'hwTkQuality': (10, -1, 9),
    'pdgId': (500, -250, 250),
    'puppiWeight': (1000, 0, 1000),
    'z0': (2000, -1, 1)
}


genmch_vars = {
    'deta': (400, -0.02, 0.02),
    'dphi': (400, -0.02, 0.02),
    'dR': (400, 0, 0.04)
}