"""
Relativistic two-body reaction kinematics
"""

import math


class TwoBody:
    def __init__(self, m1, m2, m3, m4, ek):
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4
        self.ek = ek

# defaults
        self.ncoscm = 500
        self.nogo = False

        self.cmcos3max = 2.0
        self.cmcos4max = 2.0
        self.e3atmaxang = -1.0
        self.e4atmaxang = -1.0
        self.theta3max = None
        self.theta4max = None

        self._compute()


    def _compute(self):
        # Mandelstam s
        self.s = (self.m1 + self.m2) ** 2 + 2.0 * self.m2 * self.ek
        if self.s <= 0.0:
            self.nogo = True
            return

        # initial CM momentum
        pcm2 = (
            (self.s - self.m1**2 - self.m2**2) ** 2
            - 4.0 * self.m1**2 * self.m2**2
        )
        if pcm2 < 0:
            self.nogo = True
            return

        self.pcm = math.sqrt(pcm2 / (4.0 * self.s))

        # CM rapidity
        acmratio = (
            (math.sqrt(self.m2**2 + self.pcm**2) + self.pcm)
            / self.m2
        )
        self.cmrap = math.log(acmratio)
        self.thesinh = math.sinh(self.cmrap)
        self.thecosh = math.cosh(self.cmrap)

        # final-state CM momentum
        pcmp2 = (
            (self.s - self.m3**2 - self.m4**2) ** 2
            - 4.0 * self.m3**2 * self.m4**2
        )
        if pcmp2 < 0:
            self.nogo = True
            return

        self.pcmp = math.sqrt(pcmp2 / (4.0 * self.s))

        # CM energies
        self.e03 = math.sqrt(self.pcmp**2 + self.m3**2)
        self.e04 = math.sqrt(self.pcmp**2 + self.m4**2)

        # lab-frame extrema
        self.emax3 = self.e03 * self.thecosh + self.pcmp * self.thesinh - self.m3
        self.emin3 = self.e03 * self.thecosh - self.pcmp * self.thesinh - self.m3
        self.emax4 = self.e04 * self.thecosh + self.pcmp * self.thesinh - self.m4
        self.emin4 = self.e04 * self.thecosh - self.pcmp * self.thesinh - self.m4

        # max ejectile angle
        if self.m3 > 0.0:
            thetatest = self.pcmp / (self.m3 * self.thesinh)
            
            if thetatest < 1.0:
                self.theta3max = math.asin(thetatest)

                patmax = (
                    self.e03 * math.cos(self.theta3max) * self.thesinh
                ) / (1.0 + thetatest**2 * self.thesinh**2)

                eatmax = math.sqrt(patmax**2 + self.m3**2)
                self.e3atmaxang = eatmax - self.m3
                
                self.cmcos3max = (
                    (eatmax - self.e03 * self.thecosh)
                    / (self.pcmp * self.thesinh)
                )
               
                

        # elastic symmetry case
        if (self.m1 + self.m2) == (self.m3 + self.m4):
            if abs(thetatest - 1.0) < 1e-3:
                self.theta3max = math.pi / 2.0
                self.cmcos3max = -1.0
                self.e3atmaxang = (
                    self.e03 * self.thecosh
                    + self.cmcos3max * self.pcmp * self.thesinh
                    - self.m3
                )

        # max recoil angle
        if self.m4 > 0.0:
            thetatest = self.pcmp / (self.m4 * self.thesinh)
            if thetatest < 1.0:
                self.theta4max = math.asin(thetatest)
                patmax = (
                    self.e04 * math.cos(self.theta4max) * self.thesinh
                ) / (1.0 + thetatest**2 * self.thesinh**2)
                eatmax = math.sqrt(patmax**2 + self.m4**2)
                self.e4atmaxang = eatmax - self.m4
                self.cmcos4max = (
                    (eatmax - self.e04 * self.thecosh)
                    / (self.pcmp * self.thesinh)
                )

        if (self.m1 + self.m2) == (self.m3 + self.m4):
            if abs(thetatest - 1.0) < 1e-3:
                self.theta4max = math.pi / 2.0
                self.cmcos4max = 1.0
                self.e4atmaxang = (
                    self.e04 * self.thecosh
                    - self.cmcos4max * self.pcmp * self.thesinh
                    - self.m4)
def get_points(self, kx, ky):
        """
        Generate kinematic curves.
        Returns list: [x0, y0, x1, y1, ...]
        """
        pts = []

        for i in range(-self.ncoscm, self.ncoscm + 1):
            coscm = i / self.ncoscm
            sincm = math.sqrt(max(0.0, 1.0 - coscm**2))

            ppar3 = self.pcmp * self.thecosh * coscm + self.e03 * self.thesinh
            pperp3 = self.pcmp * sincm
            ptot3 = math.hypot(ppar3, pperp3)

            ppar4 = -self.pcmp * self.thecosh * coscm + self.e04 * self.thesinh
            pperp4 = self.pcmp * sincm
            ptot4 = math.hypot(ppar4, pperp4)

            q3 = math.acos(ppar3 / ptot3) if ptot3 > 0 else 0.0
            q4 = math.acos(ppar4 / ptot4) if ptot4 > 0 else 0.0
            qcm = math.acos(coscm)

            e3 = self.e03 * self.thecosh + coscm * self.pcmp * self.thesinh - self.m3
            e4 = self.e04 * self.thecosh - coscm * self.pcmp * self.thesinh - self.m4

            v3 = ptot3 / (e3 + self.m3)
            v4 = ptot4 / (e4 + self.m4)

            x = self._select(kx, q3, q4, qcm, coscm, e3, e4, v3, v4)
            y = self._select(ky, q3, q4, qcm, coscm, e3, e4, v3, v4)

            pts.extend((x, y))

        return pts

def _select(self, k, q3, q4, qcm, coscm, e3, e4, v3, v4):
        if k == 1:
            return q3
        if k == 2:
            return q4
        if k == 3:
            return qcm
        if k == 4:
            return coscm
        if k == 5:
            return e3
        if k == 6:
            return e4
        if k == 7:
            return v3
        if k == 8:
            return v4
        return 0.0