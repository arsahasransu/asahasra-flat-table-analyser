linkvartohist = {
    'n': (51, -1, 50, 'multiplicity'),

    'mass': (500, 0, 0.5, 'mass [GeV]'),

    'e': (101, -1, 100, 'energy [GeV]'),
    'pt': (101, -1, 100, 'p_{T} [GeV]'),  # hw precision = 0.25 GeV
    'tkPt': (101, -1, 100, 'track p_{T} [GeV]'),

    'eta': (62, -3.1, 3.1, '#eta'),  # hw precision = 0.00436 GeV
    'caloeta': (62, -3.1, 3.1, 'calo. #eta'),
    'caloEta': (62, -3.1, 3.1, 'calo. #eta'),
    'tkEta': (62, -3.1, 3.1, 'track #eta'),

    'phi': (66, -3.3, 3.3, '#phi'),  # hw precision = 0.00436 GeV
    'calophi': (66, -3.3, 3.3, 'calo. #phi'),
    'caloPhi': (66, -3.3, 3.3, 'calo. #phi'),
    'tkPhi': (66, -3.3, 3.3, 'track #phi'),

    'charge': (5, -2, 3, 'charge'),
    'vz': (50, -25, 25, 'v_{z} [cm]'),
    # 'vz': (100, -100, 100, 'v_{z} [cm]'),
    'z0': (500, -100, 100, 'z0'),

    'prompt': (5, -1, 4, 'prompt status'),
    'mask': (5, -1, 4, 'boolean mask'),

    'hwQual': (10, -1, 9, 'HW qual.'),
    'hwTkQuality': (10, -1, 9, 'HW track quality'),

    'pfIso': (100000, 0, 100, 'PF rel. iso.'),
    'puppiIso': (100000, 0, 100, 'puppi rel. iso.'),
    'absPuppiIso': (100000, 0, 1000, 'puppi absolute. iso.'),
    'tkIso': (100000, 0, 100, 'track rel. iso.'),

    'pdgId': (500, -250, 250, 'pid'),
    'puppiWeight': (110, 0.5, 1.05, 'puppi weight'),

    # 'deta': (140, -7, 7, '#Delta#eta'),
    # 'dphi': (280, -7, 7, '#Delta#phi'),
    # 'dR': (280, 0, 7, '#DeltaR'),
    'deta': (200, -0.2, 0.2, '#Delta#eta'),
    'dphi': (200, -0.4, 0.4, '#Delta#phi'),
    'dR': (100, 0.0, 0.2, '#DeltaR'),

    'bin2dR': (5000, 0.0, 0.5, '#DeltaR')
}
