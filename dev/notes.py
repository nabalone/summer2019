[[0.33451957 0.06761566 0.09252669 0.50533808 0.        ]
 [0.06666667 0.33333333 0.26666667 0.33333333 0.        ]
 [0.23684211 0.23684211 0.22368421 0.30263158 0.        ]
 [0.33333333 0.05555556 0.16666667 0.44444444 0.        ]
 [0.33333333 0.         0.         0.         0.66666667]]
bal_score:0.40 diag:6
rf_1rands_50ests_a
[0.00453032 0.04792162 0.01190929 0.00330492 0.00799752 0.04186488
 0.06321143 0.05479711 0.01778947 0.02678721 0.00770866 0.01212733
 0.04920281 0.02727356 0.02216323 0.01336321 0.02533585 0.00722506
 0.02320961 0.01157646 0.00571681 0.00622712 0.01379659 0.03008411
 0.07041805 0.01651163 0.04319289 0.01629982 0.02249274 0.04295345
 0.01200805 0.01388703 0.02357222 0.01329819 0.0221316  0.02743934
 0.02231006 0.00897926 0.01169115 0.01232131 0.02070725 0.03492475
 0.00984172 0.01989531]
47.43408536911011


sort by photometric data points, top 50 or so, 
SLSN-1
sn.space


preprocessing


['KronRad (kpc)_3', 'separation (kpc)_3', 'area (kpc^2)_3', 'sep/sqrt(area) (kpc)_3', \
 'KronMag_3', 'Abs. Mag_3','Ellipticity_3', 'Z_3', 'pixelRank_3', 'chance coincidence_3',\
 'KronRad (kpc)_4', 'separation (kpc)_4', 'area (kpc^2)_4', 'sep/sqrt(area) (kpc)_4', 'x_4',\
 'y_4', 'KronMag_4', 'Abs. Mag_4', 'Angle_4', 'Ellipticity_4', 'RA_4', 'Host RA_4', 'DEC_4',\
 'Host Dec_4', 'Discrepency (arcsecs)_4', 'Z_4', 'SDSS Photoz_4', 'pixelRank_4',\
 'chance coincidence_4', 'KronRad (kpc)_5', 'separation (kpc)_5', 'area (kpc^2)_5',\
 'sep/sqrt(area) (kpc)_5', 'x_5', 'y_5', 'KronMag_5', 'Abs. Mag_5', 'Angle_5', 'Ellipticity_5',\
 'RA_5', 'Host RA_5', 'DEC_5', 'Host Dec_5', 'Discrepency (arcsecs)_5', 'Z_5', 'SDSS Photoz_5',\
 'pixelRank_5', 'chance coincidence_5', 'KronRad (kpc)_6', 'separation (kpc)_6',\
 'area (kpc^2)_6', 'sep/sqrt(area) (kpc)_6', 'x_6', 'y_6', 'KronMag_6', 'Abs. Mag_6',\
 'Angle_6', 'Ellipticity_6', 'RA_6', 'Host RA_6', 'DEC_6', 'Host Dec_6',\
 'Discrepency (arcsecs)_6', 'Z_6', 'SDSS Photoz_6', 'pixelRank_6', 'chance coincidence_6']
 
 randomly select with over 30 data points, bs, cs, and b/cs for bc. 
 get more bcs and sls, otherwise old coords
 not all will be in 3pi
 
 
 
 ../vvillar/ps1_host/fits/
 
 batch size smaller

switch to a gpu to go faster
reduce number of filters
make sure pooling after.
can batch size the evaluating
