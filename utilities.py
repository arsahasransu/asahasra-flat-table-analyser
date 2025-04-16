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
    'pt': (101, -1, 100),
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


puppiel_vars = charged_threep_vars | {
    'hwEmID': (10, -1, 9),
    'clEmEt': (101, -1, 100),
    'clPt': (101, -1, 100),
    'pfEmID': (200, -1, 1),
    'pfPuID': (200, -1, 1)
}


genmch_vars = {
    'deta': (400, -0.02, 0.02),
    'dphi': (100, 0, 0.02),
    'genqdphi': (200, -0.02, 0.02),
    'pupelqdphi': (200, -0.02, 0.02),
    'dR': (40000, 0, 4)
}