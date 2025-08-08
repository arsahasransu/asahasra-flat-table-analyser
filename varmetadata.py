linkvartohist = {
    'n': (51, -1, 50, 'multiplicity'),

    'mass': (500, 0, 0.5, 'mass [GeV]'),

    'e': (101, -1, 100, 'energy [GeV]'),
    'pt': (101, -1, 100, 'p_{T} [GeV]'),  # hw precision = 0.25 GeV
    'tkPt': (101, -1, 100, 'track p_{T} [GeV]'),

    'eta': (54, -2.7, 2.7, '#eta'),  # hw precision = 0.00436 GeV
    'caloeta': (54, -2.7, 2.7, 'calo. #eta'),
    'caloEta': (54, -2.7, 2.7, 'calo. #eta'),
    'tkEta': (54, -2.7, 2.7, 'track #eta'),

    'phi': (66, -3.3, 3.3, '#phi'),  # hw precision = 0.00436 GeV
    'calophi': (66, -3.3, 3.3, 'calo. #phi'),
    'caloPhi': (66, -3.3, 3.3, 'calo. #phi'),
    'tkPhi': (66, -3.3, 3.3, 'track #phi'),

    'charge': (5, -2, 3, 'charge'),
    'vz': (500, -100, 100, 'v_{z} [cm]'),
    'z0': (500, -100, 100, 'z0'),

    'prompt': (5, -1, 4, 'prompt status'),

    'hwQual': (10, -1, 9, 'HW qual.'),
    'hwTkQuality': (10, -1, 9, 'HQ track quality'),

    'pfIso': (100000, 0, 100, 'PF rel. iso.'),
    'puppiIso': (100000, 0, 100, 'puppi rel. iso.'),
    'tkIso': (100000, 0, 100, 'track rel. iso.'),

    'pdgId': (500, -250, 250, 'pid'),
    'puppiWeight': (110, 0, 1.1, 'puppi weight'),

    'deta': (140, -7, 7, '#Delta#eta'),
    'dphi': (140, -1, 6, '#Delta#phi'),
    'dR': (1000, 0, 10, '#DeltaR')
    # 'deta': (100, -0.2, 0.2, '#Delta#eta'),
    # 'dphi': (50, -0.1, 0.4, '#Delta#phi'),
    # 'dR': (100, 0, 0.4, '#DeltaR')
}
