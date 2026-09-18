from importlib.resources import files

import Shadow

from shadowpy.beamline import BeamLine
from shadowpy.optical_elements import Screen, Slits, ToroidalMirror
from shadowpy.sources import BendingMagnet
from shadowpy.utils import ReferenceFrame, rotation_matrix

# Beamline specification files ship as package data, so they are resolved
# relative to this package regardless of whether shadowpy is installed
# (site-packages) or run from a source checkout (editable install).
SPECS_DIR = str(files("shadowpy.beamline_config"))


# Utilize the cax-scripts repository for additional functionality
# sys.path.append(os.path.expanduser("~/repos/cax-scripts"))


class CAXSim(BeamLine):
    """
    Class to simulate the CARCARÁ-X (CAX) beamline using ShadowPy. 
    
    This class inherits from the Beamline class and provides methods to set 
    up and run simulations for the CAX beamline, perform scans and save 
    results.
    """

    CAX_SPECS = f"{SPECS_DIR}/CAX"

    sirius_frame = ReferenceFrame(name="sirius")

    R1 = rotation_matrix('x', 90)
    R2 = rotation_matrix('y', 180)
    shadow_R = R1 @ R2

    shadow_frame = sirius_frame.child_frame(relative_origin=[0, 0, 0], 
                                            relative_rotation=shadow_R, 
                                            name="shadow")

    def __init__(self, total_rays: int = 100000):
        # Initialize the elements



        self.total_rays = total_rays        
        source = BendingMagnet("B1", f"{self.CAX_SPECS}/source B1.txt")

        self.mirror = ToroidalMirror(name="M1", 
                                     specification_file=f"{self.CAX_SPECS}/mirror M1.txt")
        self.dvf_A1 = Slits(name="DVF_A1", 
                            specification_file=f"{self.CAX_SPECS}/slit_A1.txt")
        self.dvf_B1 = Screen(name="DVF_B1", 
                            specification_file=f"{self.CAX_SPECS}/dvf_B1.txt")
        # self.dvf_B1.pixel_size = 0.48/1000 # .48 μm in mm
        beam = Shadow.Beam()
        
        source.shadow_oe.NPOINT = self.total_rays

        super().__init__(name="CAX", beam=beam, source=source, 
                 optical_elements=[self.mirror, self.dvf_A1, self.dvf_B1], 
                 beamline_frame=self.shadow_frame, 
                 lab_frame=self.sirius_frame)
        
        self.trace()