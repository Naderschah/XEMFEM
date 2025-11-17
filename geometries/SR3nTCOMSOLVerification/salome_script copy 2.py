from salome.shaper import model
from SketchAPI import SketchAPI_Circle, SketchAPI_Point, SketchAPI_Arc, SketchAPI_Line
import ModelAPI
from salome.shaper import geom
import uuid
"""
Direct transcription of the COMSOL matlab code 
"""

MODELDO = True

## ------------   Constants --------------------------
from dataclasses import dataclass
import math

@dataclass
class GeometryParams:
  # Diameter of the field shaping rings
  WireDiameter: float = 0.002  # 2 [mm]
  # Radial position of the center of the field shaping rings
  WireRadialPosition: float = 0.668  # 668 [mm]
  # Vertical pitch of the field shaping rings in warm condition
  WireVerticalPitch: float = 0.022  # 22 [mm]
  # Vertical position of the top most field shaping ring in warm condition
  WireTopVerticalPosition: float = -0.023  # -23 [mm]	# WARNING: multiple values: ['-12 [mm]', '-23 [mm]', '12 [mm]', '17 [mm]']
  # Position of the top face of the gate electrode (z=0 by definition)
  GateVerticalPosition: float = 0  # 0 [mm]
  # Position of the bottom face of the anode electrode
  AnodeVerticalPosition: float = 0.008  # 8 [mm]	# WARNING: multiple values: ['26 [mm]', '8 [mm]']
  # Position of the top face of the top screen electrode
  TopScreenVerticalPosition: float = 0.055  # 55 [mm]	# WARNING: multiple values: ['51 [mm]', '55 [mm]', '56 [mm]']
  # Height of the top screen electrode frame
  TopScreenHeight: float = 0.015  # 15 [mm]
  # Width of the top screen electrode frame
  TopScreenWidth: float = 0.031  # 31 [mm]
  # Deburring radius of the edges of the top screen electrode radius
  TopScreenDeburringRadius: float = 0.0016  # 1.6 [mm]
  # Height of the anode electrode frame
  AnodeHeight: float = 0.024  # 24 [mm]	# WARNING: multiple values: ['18 [mm]', '24 [mm]']
  # Width of the anode electrode frame
  AnodeWidth: float = 0.031  # 31 [mm]
  # Deburring radius of the anode electrode frame
  AnodeDeburringRadius: float = 0.0016  # 1.6 [mm]
  # Height of the gate electrode frame
  GateHeight: float = 0.02  # 20 [mm]
  # Height of the cut-out of the gate electrode frame
  GateCutOutHeight: float = 0.011  # 11 [mm]
  # Width of the gate electrode frame
  GateWidth: float = 0.031  # 31 [mm]
  # Width of the cut-out of the gate electrode frame
  GateCutOutWidth: float = 0.01  # 10 [mm]
  # Deburring radius of the gate electrode frame
  GateDeburringRadius: float = 0.0016  # 1.6 [mm]
  # Radial position of the inner face of the anode electrode frame
  AnodeRadialPosition: float = 0.667  # 667 [mm]
  # Radial position of the inner face of the gate electrode frame
  GateRadialPosition: float = 0.667  # 667 [mm]
  # Radial position of the inner face of the top screen electrode frame
  TopScreenRadialPosition: float = 0.667  # 667 [mm]
  # Radial position of the inner face of the electrodes of the top stack
  TopStackRadialPosition: float = 0.667  # 667 [mm]
  # Spacing between wires for the top stack electriodes
  TopStackWireSpacing: float = 0.005  # 5 [mm]
  # Diameter of the wires of the top screen electrode
  TopScreenWireDiameter: float = 0.000216  # 216 [um]
  # Diameter of the wires of the gate electrode
  GateWireDiameter: float = 0.000216  # 216 [um]
  # Diameter of the wires of the electrodes of the top stack
  TopStackWireDiameter: float = 0.000216  # 216 [um]
  # Vertical position of the top face of the insulating frame of the top stack electrodes
  TopStackInsulationVerticalPosition: float = 0.055  # 55 [mm]
  # Radial position of the inner face of the insulating frame of the top stack electrodes
  TopStackInsulationRadialPosition: float = 0.6643  # 664.3 [mm]
  # Width of the insulating frame of the top stack electrodes
  TopStackInsulationWidth: float = 0.0362  # 36.2 [mm]
  # Height of the insulating frame of the top stack electrodes
  TopStackInsulationHeight: float = 0.075  # 75 [mm]
  # Distance between the end of the insulating frame and the electrode wires for the top stack electrodes (only top screen and anode)
  TopStackInsulationWireGap: float = 0.0007  # 0.7 [mm]
  # Height of the reflector of the top screen electrode, meaning the inner part of the insulating frame
  TopScreenReflectorHeight: float = 0.0332  # 33.2 [mm]
  # Height of the insulating frame of the top screen electrode
  TopScreenInsulationHeight: float = 0.0183  # 18.3 [mm]
  # Distance between the end of the insulating frame and the electrode wires for the top stack electrodes (only top screen and anode, insulation on top, electrode on bottom)
  TopStackInsulationWireGapTop: float = 0.0007  # 0.7 [mm]
  # Distance between the end of the insulating frame and the electrode wires for the top stack electrodes (only top screen and anode, insulation on bottom, electrode on top)
  TopStackInsulationWireGapBottom: float = 0.0005  # 0.5 [mm]
  # Dimension A of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionA: float = 0.0308  # 30.8 [mm]	# WARNING: multiple values: ['26.8 [mm]', '30.8 [mm]']
  # Dimension B of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionB: float = 0.0005  # 0.5 [mm]
  # Dimension C of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionC: float = 0.012  # 12 [mm]	# WARNING: multiple values: ['0.5 [mm]', '12 [mm]']
  # Dimension D of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionD: float = 0.0202  # 20.2 [mm]	# WARNING: multiple values: ['20.2 [mm]', '26.8 [mm]']
  # Dimension E of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionE: float = 0.004  # 4 [mm]	# WARNING: multiple values: ['26.8 [mm]', '4 [mm]']
  # Dimension F of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionF: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['2 [mm]', '26.8 [mm]']
  # Dimension G of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionG: float = 0.0315  # 31.5 [mm]	# WARNING: multiple values: ['26.8 [mm]', '27.5 [mm]', '31.5 [mm]']
  # Dimension H of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionH: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['2 [mm]', '26.8 [mm]']
  # Dimension I of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionI: float = 0.0222  # 22.2 [mm]	# WARNING: multiple values: ['22.2 [mm]', '26.8 [mm]']
  # Dimension J of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionJ: float = 0.0005  # 0.5 [mm]	# WARNING: multiple values: ['0.5 [mm]', '26.8 [mm]']
  # Dimension K of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionK: float = 0.0027  # 2.7 [mm]	# WARNING: multiple values: ['2.7 [mm]', '26.8 [mm]']
  # Vertical position of the bottom left corner of the insulating frame of the anode electrode frame
  AnodeInsulationVerticalPosition: float = 0.0087  # 8.7 [mm]
  # Dimension A of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionA: float = 0.026  # 26 [mm]	# WARNING: multiple values: ['26 [mm]', '26.8 [mm]']
  # Dimension B of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionB: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['2 [mm]', '26 [mm]']
  # Dimension C of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionC: float = 0.004  # 4 [mm]	# WARNING: multiple values: ['26 [mm]', '4 [mm]']
  # Dimension D of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionD: float = 0.0202  # 20.2 [mm]	# WARNING: multiple values: ['20.2 [mm]', '26 [mm]']
  # Dimension E of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionE: float = 0.012  # 12 [mm]	# WARNING: multiple values: ['12 [mm]', '26 [mm]']
  # Dimension F of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionF: float = 0.0005  # 0.5 [mm]	# WARNING: multiple values: ['0.5 [mm]', '26 [mm]']
  # Dimension G of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionG: float = 0.007  # 7 [mm]	# WARNING: multiple values: ['26 [mm]', '7 [mm]']
  # Dimension H of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionH: float = 0.0005  # 0.5 [mm]	# WARNING: multiple values: ['0.5 [mm]', '26 [mm]']
  # Dimension I of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionI: float = 0.0222  # 22.2 [mm]	# WARNING: multiple values: ['22.2 [mm]', '26 [mm]']
  # Dimension J of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionJ: float = 0.02  # 20 [mm]	# WARNING: multiple values: ['20 [mm]', '26 [mm]']
  # Dimension K of the gate insulating frame (see drawing GateInsulatingFrame)
  GateInsulationDimensionK: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['2 [mm]', '26 [mm]']
  # Vertical position of the bottom left corner of the insulating frame of the gate electrode frame
  GateInsulationVerticalPosition: float = 0.0005  # 0.5 [mm]	# WARNING: multiple values: ['- 9 [mm]', '0.5 [mm]', '8.7 [mm]']
  # Dimension L of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionL: float = 0.0255  # 25.5 [mm]	# WARNING: multiple values: ['19.5 [mm]', '25.2 [mm]', '25.5 [mm]', '27.5 [mm]']
  # Dimension A of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionA: float = 0.0165  # 16.5 [mm]	# WARNING: multiple values: ['16.5 [mm]', '4 [mm]']
  # Dimension B of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionB: float = 0.0005  # 0.5 [mm]	# WARNING: multiple values: ['0.5 [mm]', '16.5 [mm]']
  # Dimension C of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionC: float = 0.0035  # 3.5 [mm]	# WARNING: multiple values: ['16.5 [mm]', '3.5 [mm]']
  # Dimension D of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionD: float = 0.0335  # 33.5 [mm]	# WARNING: multiple values: ['16.5 [mm]', '33.5 [mm]']
  # Dimension E of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionE: float = 0.0149  # 14.9 [mm]	# WARNING: multiple values: ['14.9 [mm]', '16.5 [mm]']
  # Dimension F of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionF: float = 0.0027  # 2.7 [mm]	# WARNING: multiple values: ['16.5 [mm]', '2.7 [mm]']
  # Dimension G of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionG: float = 0.0332  # 33.2 [mm]	# WARNING: multiple values: ['16.5 [mm]', '33.2 [mm]']
  # Dimension H of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionH: float = 0.0093  # 9.3 [mm]	# WARNING: multiple values: ['16.5 [mm]', '9.3 [mm]']
  # Dimension I of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionI: float = 0.0222  # 22.2 [mm]	# WARNING: multiple values: ['16.5 [mm]', '22.2 [mm]']
  # Dimension J of the top screen insulating frame (see drawing TopScreenInsulatingFrame)
  TopScreenInsulationDimensionJ: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['16.5 [mm]', '2 [mm]']
  # Vertical position of the bottom left corner of the insulating frame of the top screen electrode frame
  TopScreenInsulationVerticalPosition: float = 0.0407  # 40.7 [mm]	# WARNING: multiple values: ['36.7 [mm]', '40.7 [mm]']
  # Vertical position of the top face of the reflecting panel
  PanelVerticalPosition: float = -0.0008  # -0.8 [mm]	# WARNING: multiple values: ['-0.8 [mm]', '0.8 [mm]']
  # Width of the reflecting panel
  PanelWidth: float = 0.003  # 3 [mm]
  # Radial position of the face of the reflecting panel facing the inner volume of the TPC
  PanelRadialPosition: float = 0.664  # 664 [mm]
  # Height of the reflecting panel
  PanelHeight: float = 1.5008  # 1500.8 [mm]
  # Shrinkage of the PTFE at the liquid xenon temperature
  ShrinkageFactor: float = 1.0  # original: '1-0.014'
  # Inner diameter of the inner cryostat
  CryostatDiameter: float = 1.46  # 1460 [mm]	# WARNING: multiple values: ['1460 [mm]', '730 [mm]']
  # Vertical position of the top copper ring
  CopperRingVerticalPosition: float = -0.025  # -25 [mm]
  # Height of the top copper ring
  CopperRingHeight: float = 0.01  # 10 [mm]
  # Width of the top copper ring
  CopperRingWidth: float = 0.021  # 21 [mm]	# WARNING: multiple values: ['21 [mm]', '21 [mm] - CopperRingModificationTemporary']
  # Radial position of the inner face of the top copper ring
  CopperRingRadialPosition: float = 0.679  # 679 [mm]	# WARNING: multiple values: ['679 [mm]', '679 [mm] +CopperRingModificationTemporary', '679 [mm] -CopperRingModificationTemporary']
  # Deburring radius of the top copper ring
  CopperRingDeburringRadius: float = 0.0015  # 1.5 [mm]
  # Height of the insulator of the top copper ring (placed between the gate and the copper ring)
  CopperRingInsulatorHeight: float = 0.005  # 5 [mm]
  # Width of the insulator of the top copper ring
  CopperRingInsulatorWidth: float = 0.028  # 28 [mm]
  # Radial position of the inner face of the insulator of the top copper ring
  CopperRingInsulatorRadialPosition: float = 0.679  # 679 [mm]
  # Vertical position of the top face of the insulator of the top copper ring
  CopperRingInsulatorVerticalPosition: float = -0.02  # -20 [mm]
  # Height of the insulator of the top copper ring (placed between the gate and the copper ring)
  CopperRingInsulationHeight: float = 0.005  # 5 [mm]
  # Radial position of the inner face of the insulator of the top copper ring
  CopperRingInsulationRadialPosition: float = 0.679  # 679 [mm]
  # Vertical position of the top face of the insulator of the top copper ring
  CopperRingInsulationVerticalPosition: float = -0.02  # -20 [mm]
  # Width of the insulator of the top copper ring
  CopperRingInsulationWidth: float = 0.028  # 28 [mm]
  # Spacing between wires for the bottom stack electriodes
  BottomStackWireSpacing: float = 0.0075  # 7.5 [mm]	# WARNING: multiple values: ['5 [mm]', '7.5 [mm]']
  # Height of the cathode electrode frame
  CathodeHeight: float = 0.02  # 20 [mm]
  # Width of the cathode electrode frame
  CathodeWidth: float = 0.024  # 24 [mm]
  # Deburring radius of the cathode electrode frame
  CathodeDeburringRadius: float = 0.0005  # 0.5 [mm]
  # Radius of the rounding of the cathode electrode frame
  CathodeRoundingRadius: float = 0.01  # 10 [mm]
  # Vertical position of the upper face of the cathode electrode frame
  CathodeVerticalPosition: float = -1.5028  # -1502.8 [mm]
  # Diameter of the wires of the cathode electrode
  CathodeWireDiameter: float = 0.0003  # 300 [um]
  # Radial position of the cathode electrode frame
  CathodeRadialPosition: float = 0.6735  # 673.5 [mm]
  # Radial position of the inner surface of the bottom screen electrode frame
  BottomScreenRadialPosition: float = 0.6725  # 672.5 [mm]
  # Width of the bottom screen electrode frame
  BottomScreenWidth: float = 0.025  # 25 [mm]
  # Hieght of the bottom screen electrode frame
  BottomScreenHeight: float = 0.015  # 15 [mm]
  # Deburring radius of the bottom screen electrode frame
  BottomScreenDeburringRadius: float = 0.0005  # 0.5 [mm]
  # Rounding radious of the bottom screen electrode frame
  BottomScreenRoundingRadius: float = 0.0075  # 7.5 [mm]
  # Vertical position of the lower face of the bottom screen electrode frame
  BottomScreenVerticalPosition: float = -1.558  # -1558 [mm]	# WARNING: multiple values: ['-1558 [mm]', '1558 [mm]']
  # Diameter of the wires of the bottom screen electrode
  BottomScreenWireDiameter: float = 0.000216  # 216 [um]
  # Dimension A of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionA: float = 0.01  # 10 [mm]
  # Dimension B of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionB: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['10 [mm]', '2 [mm]']
  # Dimension C of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionC: float = 0.012  # 12 [mm]	# WARNING: multiple values: ['10 [mm]', '12 [mm]']
  # Dimension D of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionD: float = 0.007  # 7 [mm]	# WARNING: multiple values: ['10 [mm]', '7 [mm]']
  # Dimension E of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionE: float = 0.008  # 8 [mm]	# WARNING: multiple values: ['10 [mm]', '8 [mm]']
  # Dimension F of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionF: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['10 [mm]', '2 [mm]']
  # Dimension G of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionG: float = 0.012  # 12 [mm]	# WARNING: multiple values: ['10 [mm]', '12 [mm]']
  # Dimension H of the cathode insulating frame (see drawing CathodeInsulatingFrame)
  CathodeInsulationDimensionH: float = 0.005  # 5 [mm]	# WARNING: multiple values: ['12 [mm]', '5 [mm]']
  # Radial position of the insulating frame of the cathode electrode frame
  CathodeInsulationRadialPosition: float = 0.6855  # 685.5 [mm]
  # Vertical position of the insulating frame of the cathode electrode frame
  CathodeInsulationVerticalPosition: float = -1.5028  # -1502.8 [mm]	# WARNING: multiple values: ['-1502.8 [mm]', '1502.8 [mm]']
  # Height of the bottom stack reflector
  BottomStackInsulatorHeight: float = 0.05378  # 53.78 [mm]
  # Width of the bottom stack reflector+
  BottomStackInsulationWidth: float = 0.0042  # 4.2 [mm]
  # Height of the bottom stack reflector
  BottomStackInsulationHeight: float = 0.05378  # 53.78 [mm]
  # Vertical position of the top face of the bottom stack reflector
  BottomStackInsulationVerticalPosition: float = -1.5035  # -1503.5 [mm]	# WARNING: multiple values: ['-1503.5 [mm]', '1503.5 [mm]']
  # Radial position of the inner surface of the bottom stack reflector
  BottomStackInsulationRadialPosition: float = 0.664  # 664 [mm]
  # Diameter of the grooves of the bottom stack reflector
  BottomStackInsulationGrooveDiameter: float = 0.001  # 1 [mm]	# WARNING: multiple values: ['1 [mm]', '4 [mm]']
  # Vertical position of the top most groove of the bottom stack reflector with respect to its upper surface
  BottomStackInsulationGrooveTopVerticalPosition: float = 0.00178  # 1.78 [mm]
  # Vertical pitch between two grooves in the botto stack reflector+
  BottomStackInsualtionGrooveVerticalPitch: float = 0.002  # 2 [mm]
  # Vertical pitch between two grooves in the botto stack reflector+
  BottomStackInsulationGrooveVerticalPitch: float = 0.002  # 2 [mm]
  # Height of the field shaping guard
  GuardHeight: float = 0.015  # 15 [mm]
  # Width of the field shaping guard
  GuardWidth: float = 0.005  # 5 [mm]	# WARNING: multiple values: ['10 [mm]', '5 [mm]']
  # Vertical position of the center of the top most field shaping guard
  GuardTopVerticalPosition: float = -0.078125  # -78.125 [mm]	# WARNING: multiple values: ['-63.125 [mm]', '-78.125 [mm]']
  # Radial position of the inner face of the field shaping guard
  GuardRadialPosition: float = 0.677687  # 677.687 [mm]
  # Vertical pitch of the field shaping guard
  GuardVerticalPitch: float = 0.022  # 22 [mm]
  # Number of field shaping guards
  GuardNumber: int = 64  # original: '64'
  # Radial position of the bottom holder of the bottom screen electrode frame
  BottomScreenInsulationHolderRadialPosition: float = 0.69025  # 690.25 [mm]
  # Radial position of the bottom holder of the bottom screen electrode frame
  BottomScreenHolderRadialPosition: float = 0.69025  # 690.25 [mm]
  # Vertical position of the bottom holder of the bottom screen electrode frame
  BottomScreenHolderVerticalPosition: float = -1.558  # -1558 [mm]	# WARNING: multiple values: ['-1558 [mm]', '690.25 [mm]']
  # Width of the bottom holder of the bottom screen electrode frame
  BottomScreenHolderWidth: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['-1558 [mm]', '2 [mm]']
  # Height of the bottom holder of the bottom screen electrode frame
  BottomScreenHolderHeight: float = 0.005  # 5 [mm]	# WARNING: multiple values: ['2 [mm]', '5 [mm]']
  # Chamfer of the bottom holder of the bottom screen electrode frame
  BottomScreenHolderChamfer: float = 0.001  # 1 [mm]	# WARNING: multiple values: ['1 [mm]', '5 [mm]']
  # Position of the closest wire to the guard frame
  GateFirstWireRadialPosition: float = 0.0  # original: 'TopStackRadialPosition-4.5[mm]'	# WARNING: multiple values: ['TopStackRadialPosition-4.497[mm]', 'TopStackRadialPosition-4.5[mm]']
  # Position of the closest wire to the anode frame
  AnodeFirstWireRadialPosition: float = 0.0  # original: 'TopStackRadialPosition-7[mm]'	# WARNING: multiple values: ['TopStackRadialPosition-7[mm]', 'TopStackRadialPosition-8.101[mm]']
  # Position of the closest wire to the top screen frame
  TopScreenFirstWireRadialPosition: float = 0.0  # original: 'TopStackRadialPosition-4.5[mm]'
  # Position of the closest wire to the cathode frame
  CathodeFirstWireRadialPosition: float = 0.6675  # 667.5 [mm]	# WARNING: multiple values: ['667.5 [mm]', 'TopStackRadialPosition-7[mm]']
  # Position of the closest wire to the frame for the bottom stack electrodes (cathode and bottom screen)
  BottomStackFirstWireRadialPosition: float = 0.6675  # 667.5 [mm]
  # Height of the bottom PMTs array reflector
  BottomPMTReflectorHeight: float = 0.008  # 8 [mm]
  # Height of the bottom PMTs array reflector
  BottomPMTsReflectorHeight: float = 0.008  # 8 [mm]
  # Diameter of the bottom PMTs array reflector
  BottomPMTsArrayReflectorDiameter: float = 1.395  # 1395 [mm]
  # Height of the bottom PMTs array reflector
  BottomPMTsArrayReflectorHeight: float = 0.008  # 8 [mm]
  # Vertical position of the top face of the bottom PMTs array reflector
  BottomPMTsArrayReflectorVerticalPosition: float = -1.56  # -1560 [mm]	# WARNING: multiple values: ['-1560 [mm]', '1560 [mm]']
  # Width of the notch in the bottom PMTs array reflector
  BottomPMTsArrayReflectorNotchWidth: float = 0.002  # 2 [mm]
  # Radial position of the notch in the bottom PMTs array reflector
  BottomPMTsArrayReflectorNotchRadialPosition: float = 0.69025  # 690.25 [mm]	# WARNING: multiple values: ['2 [mm]', '690.25 [mm]']
  # Depth of the notch in the bottom PMTs array reflector
  BottomPMTsArrayReflectorNotchDepth: float = 0.003  # 3 [mm]	# WARNING: multiple values: ['2 [mm]', '3 [mm]']
  # Diameter of the holes for the PMT in the bottom PMT array
  BottomPMTsArrayReflectorPMTHoleDiameter: float = 0.069  # 69 [mm]
  # Chamfer of the holes for the PMT in the bottom PMT array
  BottomPMTsArrayReflectorPMTHoleChamfer: float = 0.002  # 2 [mm]	# WARNING: multiple values: ['2 [mm]', '69 [mm]']
  # Depth of the holes for the PMT in the bottom PMT array
  BottomPMTsArrayReflectorPMTHoleDepth: float = 0.0029  # 2.9 [mm]	# WARNING: multiple values: ['2 [mm]', '2.9 [mm]']
  # Position of the outest PMT in the hexagonal bottom PMTs array
  BottomPMTsArrayReflectorPMTFirstRadialPosition: float = 0.648  # 648 [mm]	# WARNING: multiple values: ['2.9 [mm]', '648 [mm]']
  # Position of the outest PMT in the hexagonal bottom PMTs array
  BottomPMTsArrayPMTFirstRadialPosition: float = 0.648  # 648 [mm]
  # Diameter of the holes for the PMT in the bottom PMT array
  BottomPMTsArrayReflectorPMTHoleHousingRadius: float = 0.0068  # 6.8 [mm]	# WARNING: multiple values: ['6.8 [mm]', '69 [mm]']
  # Diameter of the hole housing the PMT in the bottom PMT array
  BottomPMTsArrayReflectorPMTHoleHousingDiameter: float = 0.0785  # 78.5 [mm]	# WARNING: multiple values: ['6.8 [mm]', '78.5 [mm]']
  # Radial pitch of the holes for the PMT in the bottom PMT array
  BottomPMTsArrayReflectorPMTHoleRadialPitch: float = 0.081  # 81 [mm]	# WARNING: multiple values: ['69 [mm]', '81 [mm]']
  # Major diameter of the PMT
  PMTMajorDiameter: float = 0.076  # 76 [mm]
  # Quarzdiameter of the PMT
  PMTQuarzDiameter: float = 0.076  # 76 [mm]
  # Quartz diameter of the PMT
  PMTQuartzDiameter: float = 0.076  # 76 [mm]
  # Diameter of the body of the PMT
  PMTBodyDiameter: float = 0.0533  # 53.3 [mm]
  # Diameter of the basis of the PMT
  PMTBasisDiameter: float = 0.0318  # 31.8 [mm]
  # Height of the quartz of the PMT
  PMTQuartzHeight: float = 0.0135  # 13.5 [mm]
  # Height of the rounding of the PMT body
  PMTRoundingHeight: float = 0.0271  # 27.1 [mm]
  # Height of the body of the PMT (including rounding quartz)
  PMTBodyHeight: float = 0.1005  # 100.5 [mm]	# WARNING: multiple values: ['100.5[mm]', '114[mm]', '73.4 [mm]']
  # Height of the PMT basis
  PMTBasisHeight: float = 0.00955  # 9.55 [mm]
  # Position of the outest PMT in the hexagonal top PMTs array
  TopPMTsArrayPMTFirstRadialPosition: float = 0.648  # 648 [mm]
  # Diameter of the top PMTs array reflector
  TopPMTsArrayReflectorDiameter: float = 1.412  # 1412 [mm]
  # Height of the top PMTs array reflector
  TopPMTsArrayReflectorHeight: float = 0.008  # 8 [mm]
  # Depth of the notch in the bottom PMTs array reflector
  TPMTsArrayReflectorNotchDepth: float = 0.003  # 3 [mm]
  # Chamfer of the holes for the PMT in the top PMT array
  TopPMTsArrayReflectorPMTHoleChamfer: float = 0.002  # 2 [mm]
  # Depth of the holes for the PMT in the top PMT array
  TopPMTsArrayReflectorPMTHoleDepth: float = 0.0029  # 2.9 [mm]
  # Diameter of the holes for the PMT in the top PMT array
  TopPMTsArrayReflectorPMTHoleDiameter: float = 0.069  # 69 [mm]
  # Diameter of the hole housing the PMT in the top PMT array
  TopPMTsArrayReflectorPMTHoleHousingDiameter: float = 0.0785  # 78.5 [mm]
  # Radial pitch of the holes for the PMT in the top PMT array
  TopPMTsArrayReflectorPMTHoleRadialPitch: float = 0.081  # 81 [mm]
  # Vertical position of the top face of the top PMTs array reflector
  TopPMTsArrayReflectorVerticalPosition: float = 0.0739  # 73.9 [mm]	# WARNING: multiple values: ['69.9 [mm]', '69.96 [mm]', '73.9 [mm]']
  # Outer diameter of the diving bell
  BellOuterDiameter: float = 0.713  # 713 [mm]
  # Height of the diving bell
  BellHeight: float = 0.264  # 264 [mm]
  # Minor thickness of the diving bell
  BellMinorThickness: float = 0.004  # 4 [mm]
  # Major thickness of the diving bell
  BellMajorThickness: float = 0.005  # 5 [mm]
  # Vertical position of the increase in thickness of the diving bell
  BellThicknessIncreaseVerticalPosition: float = 0.159  # 159 [mm]
  # Vertical position of the top face of the bell
  BellVerticalPosition: float = 0.239  # 239 [mm]
  # Vertical position of the bottom of the straight section of the inner cryostat
  CryostatStraightVerticalPosition: float = -1.60824  # -1608.24 [mm]
  # Height of the straight section of the inner cryostat
  CryostatStraightHeight: float = 1.8671  # 1867.1 [mm]
  # Liquid level above the gate
  LiquidLevel: float = 0.004  # 4 [mm]
  # LowerCryostataUpperRadius (no description found)
  LowerCryostataUpperRadius: float = 0.22  # 220 [mm]
  # Radius of the upper part of the bottom of the cryostat
  LowerCryostatUpperRadius: float = 0.22  # 220 [mm]
  # Radius of the lower part of the bottom of the cryostat
  CryostatBottomLowerRadius: float = 1.2  # 1200 [mm]	# WARNING: multiple values: ['1200 [mm]', '220 [mm]']
  # Radius of the upper part of the bottom of the cryostat
  CryostatBottomUpperRadius: float = 0.22  # 220 [mm]
  # Vertical position of the center of the lower roudning of the bottom of the inner cryostat
  CryostatBottomLowerCenterVerticalPosition: float = -0.760142  # -760.142 [mm]
  # Radius of the lower part of the top of the cryostat
  CryostatTopLowerRadius: float = 0.22  # 220 [mm]	# WARNING: multiple values: ['1200 [mm]', '220 [mm]']
  # Radius of the upper part of the top of the cryostat
  CryostatTopUpperRadius: float = 1.2  # 1200 [mm]	# WARNING: multiple values: ['1200 [mm]', '220 [mm]']
  # Diameter of the pipe opening on the top of the inner cryostat
  CryostatTopPipeDiameter: float = 0.378  # 378 [mm]
  # dV (no description found)
  dV: float = 0.0
  # VoltageDropFieldCage (no description found)
  VoltageDropFieldCage: float = 0.0
  # Vertical position of the top wire of the field cage in warm condition
  WireTopFieldCageVerticalPosition: float = 0.0  # original: 'WireTopVerticalPosition-2*WireVerticalPitch'
  # TopExtraction (no description found)
  TopExtraction: float = 0.0
  # TopExport (no description found)
  TopExport: float = 0.0
  # BottomExtraction (no description found)
  BottomExtraction: float = 0.0
  # BottomExport (no description found)
  BottomExport: float = 0.0
  # RadialExtraction (no description found)
  RadialExtraction: float = 0.0
  # RadialExport (no description found)
  RadialExport: float = 0.0
  # Change of the radial width and position of the copper ring - TEMPORARY
  CopperRingModificationTemporary: float = 0.002  # 2 [mm]
  # Dimension L of the anode insulating frame (see drawing AnodeInsulatingFrame)
  AnodeInsulationDimensionN: float = 0.0148  # 14.8 [mm]	# WARNING: multiple values: ['13 [mm]', '14.8 [mm]', '19.5 [mm]']
  # PMTCasing (no description found)
  PMTCasing: float = 0.0
  # bPMTCasing (no description found)
  bPMTCasing: float = 0.0
  # ResistanceForFirstAndLastSetofDividers (no description found)
  ResistanceForFirstAndLastSetofDividers: float = 0.0
  # RFirst5FieldShaping (no description found)
  RFirst5FieldShaping: float = 0.0
  # ResistanceLastFieldShapingAndCathode (no description found)
  ResistanceLastFieldShapingAndCathode: float = 0.0
  # RInParallelFieldShaping (no description found)
  RInParallelFieldShaping: float = 0.0
  # ResistanceCircuitSplitter (no description found)
  ResistanceCircuitSplitter: float = 0.0
  # RLastTwoFieldShaping (no description found)
  RLastTwoFieldShaping: float = 0.0
  # RFieldShapingToCathode (no description found)
  RFieldShapingToCathode: float = 0.0
  # CopperRing (no description found)
  CopperRing: float = 0.0
  # VCopperRing (no description found)
  VCopperRing: float = 0.0
  # TopFieldShapingRing (no description found)
  TopFieldShapingRing: float = 0.0
  # VTopFieldShapingRing (no description found)
  VTopFieldShapingRing: float = 0.0
  # GateElectrode (no description found)
  GateElectrode: float = 0.0
  # VGateElectrode (no description found)
  VGateElectrode: float = 0.0
  # AnodeElectrode (no description found)
  AnodeElectrode: float = 0.0
  # VAnodeElectrode (no description found)
  VAnodeElectrode: float = 0.0
  # TopShieldingElectrode (no description found)
  TopShieldingElectrode: float = 0.0
  # VTopShieldingElectrode (no description found)
  VTopShieldingElectrode: float = 0.0
  # VPMTCasing (no description found)
  VPMTCasing: float = 0.0
  # Bell (no description found)
  Bell: float = 0.0
  # VBell (no description found)
  VBell: float = 0.0
  # ParallelSegmentGuard (no description found)
  ParallelSegmentGuard: float = 0.0
  # RParallelSegmentGuard (no description found)
  RParallelSegmentGuard: float = 0.0
  # ParallelSegmentRings (no description found)
  ParallelSegmentRings: float = 0.0
  # RParallelSegmentRings (no description found)
  RParallelSegmentRings: float = 0.0
  # CathodeToGate (no description found)
  CathodeToGate: float = 0.0
  # CathodeToGateDistance (no description found)
  CathodeToGateDistance: float = 0.0
  # V2ndFieldShaping (no description found)
  V2ndFieldShaping: float = 0.0
  # VDropFirst5FieldShaping (no description found)
  VDropFirst5FieldShaping: float = 0.0



def build_lower_cryostat(part_doc, p, makeface):
    """
    2D cross-section of the lower cryostat as ONE closed contour:
      - 3 lines: axis, top, right wall
      - lower rounding arc (arc11)
      - upper rounding arc (arc10) connecting arc11 to the right wall
      - closing line from axis down to arc11
    """

    Sketch_1 = model.addSketch(part_doc, model.defaultPlane("XOY"))
    Sketch_1.setName("LowerCryostatSketch")


    # --------------------------------------------------------
    # Levels
    # --------------------------------------------------------
    y1 = p.CryostatStraightVerticalPosition
    y2 = p.CryostatStraightVerticalPosition + abs(p.CryostatStraightVerticalPosition) + p.LiquidLevel
    R = p.CryostatDiameter / 2.0
    # --------------------------------------------------------
    # Straight section: 3 lines (no rectangle primitive)
    #   - axis:   (0, y1) -> (0, y2)
    #   - top:    (0, y2) -> (R, y2)
    #   - right:  (R, y2) -> (R, y1)
    # --------------------------------------------------------
    line_axis_rect = Sketch_1.addLine(0.0, y1, 0.0, y2)
    line_axis_rect.setName("AxisRect")
    line_top = Sketch_1.addLine(0.0, y2, R,   y2)
    line_top.setName("TopLine")
    line_right = Sketch_1.addLine(R,   y2, R, y1)
    line_right.setName("RightWall")

    # --------------------------------------------------------
    # LOWER rounding: arc11 (as in your original code)
    # --------------------------------------------------------
    R_upper = p.CryostatBottomUpperRadius
    R_lower = p.CryostatBottomLowerRadius

    term_radius_diff = R_lower - R_upper
    dx = R - R_upper
    alpha = math.asin(dx / term_radius_diff)      # sector angle [rad]

    term1 = (R_lower - R_upper) ** 2
    term2 = (R - R_upper) ** 2
    dy = math.sqrt(term1 - term2)

    cx11 = 0.0
    cy11 = p.CryostatStraightVerticalPosition + dy
    r11  = R_lower

    theta_start = -math.pi / 2.0
    theta_end   = theta_start + alpha

    xs = cx11 + r11 * math.cos(theta_start)
    ys = cy11 + r11 * math.sin(theta_start)
    xe = cx11 + r11 * math.cos(theta_end)
    ye = cy11 + r11 * math.sin(theta_end)

    arc11 = Sketch_1.addArc(cx11, cy11, xs, ys, xe, ye, False)
    arc11.setName("LowerCryostatArc")

    # --------------------------------------------------------
    # Closing line on axis: from straight section down to arc11
    #   (0, y1) -> (xs, ys)
    # --------------------------------------------------------
    line_axis = Sketch_1.addLine(0.0, y1, xs, ys)
    line_axis.setName("AxisClosingLine")

    # --------------------------------------------------------
    # UPPER rounding: replace full circle with an arc
    #  - same radius R_upper
    #  - connect arc11 end to right wall at (R, y1)
    # --------------------------------------------------------
    cx10 = R - R_upper
    cy10 = y1
    r10  = R_upper

    # Right-wall / rounding junction
    xR = R
    yR = y1

    # Upper rounding arc: from arc11 end (xe, ye) to right-wall foot (R, y1)
    # Note: we rely on geometry choice so that (xe, ye) lies on this circle.
    arc10 = Sketch_1.addArc(cx10, cy10, xe, ye, xR, yR, False)
    arc10.setName("UpperRoundingArc")

    # --------------------------------------------------------
    # Coincidence constraints to ensure ONE closed wire
    #  1) axis-rectangle to top line
    #  2) top line to right wall
    #  3) right wall to upper arc
    #  4) upper arc to lower arc
    #  5) lower arc to axis closing line
    #  6) axis closing line to axis-rectangle
    # --------------------------------------------------------
    # 1) top-left corner: AxisRect top = TopLine start
    Sketch_1.setCoincident(line_axis_rect.endPoint(), line_top.startPoint())

    # 2) top-right corner: TopLine end = RightWall start
    Sketch_1.setCoincident(line_top.endPoint(), line_right.startPoint())

    # 3) bottom-right corner: RightWall end = UpperRoundingArc end (R, y1)
    Sketch_1.setCoincident(line_right.endPoint(), arc10.endPoint())

    # 4) junction between upper and lower rounding: arc10 start = arc11 end
    Sketch_1.setCoincident(arc10.startPoint(), arc11.endPoint())

    # 5) junction between lower arc and axis-closing line: arc11 start = AxisClosingLine end
    Sketch_1.setCoincident(arc11.startPoint(), line_axis.endPoint())

    # 6) junction between axis-closing line and axis-rectangle: AxisClosingLine start = AxisRect start
    Sketch_1.setCoincident(line_axis.startPoint(), line_axis_rect.startPoint())

    if MODELDO: model.do() # finalize the sketch so its wires/faces exist
    if makeface:
      LXeCryostat = model.addFace(part_doc, [Sketch_1.result()])
      LXeCryostat.setName("LXeCryostat")
      return LXeCryostat
    else:
      return Sketch_1


from salome.shaper import model
import math

def build_upper_cryostat(part_doc, p, makeface, Sketch=None):
    """
    Upper cryostat + straight section as ONE closed contour, following the
    same construction logic as the COMSOL c12/c13 definitions.
    """

    Sketch_1 = model.addSketch(part_doc, model.defaultPlane("XOY"))
    Sketch_1.setName("UpperCryostatSketch")

    R = p.CryostatDiameter / 2.0

    # ------------------------------------------------------------------
    # Straight cylindrical section
    # ------------------------------------------------------------------
    y1 = p.LiquidLevel
    y2 = p.CryostatStraightHeight + p.CryostatStraightVerticalPosition

    # axis from liquid level y1 to cap base y2
    line_axis_rect = Sketch_1.addLine(0.0, y1, 0.0, y2)

    # bottom flat AT LIQUID LEVEL
    line_top = Sketch_1.addLine(0.0, y1, R, y1)

    # right wall from liquid level up to cap base
    line_right = Sketch_1.addLine(R, y1, R, y2)
    # ------------------------------------------------------------------
    # Parameters from COMSOL-style construction
    #   - R_lower_top:  radius of circle on the outer side (c12)
    #   - R_upper_top:  radius of circle on the axis (c13)
    # ------------------------------------------------------------------
    R_lower_top = p.CryostatTopLowerRadius   # COMSOL c12 radius
    R_upper_top = p.CryostatTopUpperRadius   # COMSOL c13 radius

    dx     = R - R_lower_top
    deltaR = R_upper_top - R_lower_top

    # COMSOL angle:
    alpha = math.asin(dx / deltaR)
    # vertical offset between c13 center and y2:
    dy    = math.sqrt(deltaR**2 - dx**2)

    # ------------------------------------------------------------------
    # Circle c13 (on axis): center (0, y2 - dy), radius = R_upper_top
    # This corresponds to the COMSOL final 'pos' for c13.
    # We only use a sector of this circle as an arc.
    # ------------------------------------------------------------------
    cx13 = 0.0
    cy13 = y2 - dy
    r13  = R_upper_top

    # COMSOL:
    #   angle = alpha
    #   rot   = 90° - alpha
    # → in radians:
    theta_s_small = math.pi/2.0 - alpha   # "rot"
    theta_e_small = math.pi/2.0           # rot + angle

    # Endpoints of the small arc on c13
    x_small_start = cx13 + r13 * math.cos(theta_s_small)  # common point with c12
    y_small_start = cy13 + r13 * math.sin(theta_s_small)
    x_small_end   = cx13 + r13 * math.cos(theta_e_small)  # axis-side point
    y_small_end   = cy13 + r13 * math.sin(theta_e_small)

    arc_small = Sketch_1.addArc(
        cx13, cy13,
        x_small_start, y_small_start,  # start: common point
        x_small_end,   y_small_end,    # end: axis-side point
        False
    )
    # If this still turns the wrong way, flip the last argument to True.

    # ------------------------------------------------------------------
    # Circle c12 (outer side): center (R - R_lower_top, y2),
    # radius = R_lower_top. We use only the arc between:
    #   - right-wall top corner (R, y2)
    #   - SAME common point as for c13 (x_small_start, y_small_start)
    # ------------------------------------------------------------------
    cx12 = R - R_lower_top
    cy12 = y2
    r12  = R_lower_top

    xR = R
    yR = y2

    arc_big = Sketch_1.addArc(
        cx12, cy12,
        xR, yR,                     # start at right-wall corner
        x_small_start, y_small_start,  # end at common point with small arc
        False
    )
    # Again, if orientation is inverted, toggle the last boolean.

    # ------------------------------------------------------------------
    # Axis-cap line: connect axis top (0, y2) to axis-side point of small arc
    # ------------------------------------------------------------------
        # axis-cap line from axis top y2 to small-arc axis-side point
    line_axis_cap = Sketch_1.addLine(0.0, y2, x_small_end, y_small_end)

    # ------------------------------------------------------------------
    # Coincidence constraints for one closed loop:
    # ------------------------------------------------------------------
    # 1) axis top ↔ axis_cap start
    Sketch_1.setCoincident(line_axis_rect.endPoint(),  line_axis_cap.startPoint())

    # 2) axis_cap end ↔ small arc end (axis side)
    Sketch_1.setCoincident(line_axis_cap.endPoint(),   arc_small.endPoint())

    # 3) small arc start ↔ big arc end
    Sketch_1.setCoincident(arc_small.startPoint(),     arc_big.endPoint())

    # 4) big arc start ↔ right wall top (R, y2)
    Sketch_1.setCoincident(arc_big.startPoint(),       line_right.endPoint())

    # 5) right wall bottom (R, y1) ↔ bottom flat right end
    Sketch_1.setCoincident(line_right.startPoint(),    line_top.endPoint())

    # 6) bottom flat left end (0, y1) ↔ axis bottom
    Sketch_1.setCoincident(line_top.startPoint(),      line_axis_rect.startPoint())

    if MODELDO: model.do()
    if makeface:
      GXeCryostat = model.addFace(part_doc, [Sketch_1.result()])
      GXeCryostat.setName("GXeCryostat")
      return GXeCryostat
    else:
      return Sketch_1


def create_polygon_from_corners(Sketch, pts):
    """
    pts: list of N unique (x, y) points in order around the boundary.
    Returns N SketchLine objects forming a closed polygon.
    """
    n = len(pts)
    lines = []

    # Create N segments: from pts[i] to pts[(i+1) % n]
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        l = Sketch.addLine(x1, y1, x2, y2)
        lines.append(l)

    # Coincidence at each vertex
    for i in range(n):
        Sketch.setCoincident(lines[i].endPoint(),
                             lines[(i + 1) % n].startPoint())

    return lines

def make_anchors(Sketch, lines, pts, line_indices):
    """
    Anchors a geometry by attaching fixed construction points to line midpoints.
    Assumes straight lines and that each line connects two points in `pts`.

    Makes auxiliary (construction) points so they are *not* part of the solid geometry.
    """
    n_pts = len(pts)
    anchors = []
    for line_idx in line_indices:
        # indices of the endpoints of the line
        i0 = line_idx
        i1 = (line_idx + 1) % n_pts

        x_mid = 0.5 * (pts[i0][0] + pts[i1][0])
        y_mid = 0.5 * (pts[i0][1] + pts[i1][1])

        # Create anchor as construction geometry
        anchor = Sketch.addPoint(x_mid, y_mid)
        anchor.setAuxiliary(True)  # <-- Mark it as construction geometry

        # Coincident constraint pins it to the line
        Sketch.setCoincident(anchor, lines[line_idx].result())

        # Fix its position for solver stabilization
        Sketch.setFixed(anchor)

        anchors.append(anchor)

    return anchors



def build_gate(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("Gate")
    else:
      dont_return = True
    # The Gate ring 
    pts = [
        [p.TopStackRadialPosition, p.GateVerticalPosition-p.GateHeight+p.GateCutOutHeight],
        [p.TopStackRadialPosition, p.GateVerticalPosition],
        [p.TopStackRadialPosition+p.GateWidth, p.GateVerticalPosition],
        [p.TopStackRadialPosition+p.GateWidth, p.GateVerticalPosition-p.GateHeight],
        [p.TopStackRadialPosition+p.GateCutOutWidth, p.GateVerticalPosition-p.GateHeight],
        [p.TopStackRadialPosition+p.GateCutOutWidth, p.GateVerticalPosition-p.GateHeight+p.GateCutOutHeight],
    ]

    lines = create_polygon_from_corners(Sketch, pts)

    anchors = make_anchors(Sketch, lines, pts, line_indices=[0, 1])

    # Fillet all external corners: indices 0,1,2,3,4
    fillet_indices = [0, 1, 2, 3, 4]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.GateDeburringRadius)


    # Gate wires 
    for i in range(floor(p.TopStackRadialPosition/p.TopStackWireSpacing)):
      # First center coordinates then a point the circle passes
      Sketch.addCircle(p.GateFirstWireRadialPosition+p.TopStackWireSpacing/4+i*p.TopStackWireSpacing, 
                       p.GateVerticalPosition,
                       p.GateFirstWireRadialPosition+p.TopStackWireSpacing/4+i*p.TopStackWireSpacing, 
                       p.GateVerticalPosition+50e-6)

    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        Gate = model.addFace(part_doc, [Sketch.result()])
        Gate.setName("Gate")
        return Gate
      else:
        return Sketch

def build_anode(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("Anode")
    else:
      dont_return = True
    x0 = p.TopStackRadialPosition
    y0 = p.AnodeVerticalPosition
    W  = p.AnodeWidth
    H  = p.AnodeHeight

    pts = [
        (x0,     y0),       # bottom-left
        (x0,     y0 + H),   # top-left
        (x0 + W, y0 + H),   # top-right
        (x0 + W, y0),       # bottom-right
    ]

    lines = create_polygon_from_corners(Sketch, pts)
    # It wouldnt anchor with 2 for some reason
    anchors = make_anchors(Sketch, lines, pts, line_indices=[0, 1, 2, 3]) 

    # Fillet all external corners: indices 0,1,2,3,4
    fillet_indices = [0,1,2,3]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.AnodeDeburringRadius)


    # Gate wires 
    for i in range(floor(p.TopStackRadialPosition/p.TopStackWireSpacing)):
      # First center coordinates then a point the circle passes
      Sketch.addCircle(p.AnodeFirstWireRadialPosition+p.TopStackWireSpacing/4+i*p.TopStackWireSpacing, 
                       p.AnodeVerticalPosition,
                       p.AnodeFirstWireRadialPosition+p.TopStackWireSpacing/4+i*p.TopStackWireSpacing, 
                       p.AnodeVerticalPosition+p.TopStackWireDiameter/2)

    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        Anode = model.addFace(part_doc, [Sketch.result()])
        Anode.setName("Anode")
        return Anode
      else:
        return Sketch

def build_cathode(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("Cathode")
    else:
      dont_return = True
    x0 = p.CathodeRadialPosition
    y0 = (p.CathodeVerticalPosition- p.CathodeHeight)*p.ShrinkageFactor
    W  = p.CathodeWidth
    H  = p.CathodeHeight

    pts = [
        (x0,     y0),       # bottom-left
        (x0,     y0 + H),   # top-left
        (x0 + W, y0 + H),   # top-right
        (x0 + W, y0),       # bottom-right
    ]

    pts = [
      (x0,     y0),       # bottom-left
      (x0,     y0 + H),   # top-left (end og rect)
      (p.CathodeRadialPosition + 0.0045, y0 + H),# to inner edge
      (p.CathodeRadialPosition + 0.0045, y0 + H +0.003),# Up 
      (p.CathodeRadialPosition + 0.0045 + 0.002, y0 + H +0.003),# To inner edge
      (p.CathodeRadialPosition + 0.0045 + 0.002, y0 + H +0.003 +0.0018),# To upper left edge 
      (p.CathodeRadialPosition + 0.0045 + 0.002 +0.013, y0 + H +0.003 +0.0018),# To right upper edge 
      (p.CathodeRadialPosition + 0.0045 + 0.015, y0 + H),# To right inner edge 
      (x0 + W, y0 + H),   # top-right (og rect)
      (x0 + W, y0),       # bottom-right

    ]

    lines = create_polygon_from_corners(Sketch, pts)
    # It wouldnt anchor with 2 for some reason
    anchors = make_anchors(Sketch, lines, pts, line_indices=[i for i in range(len(lines))]) 

    # Fillet all external corners: indices 0,1,2,3,4
    fillet_indices = [1,-2]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.CathodeDeburringRadius)

    fillet_indices = [0,-1]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.CathodeRoundingRadius)

    fillet_indices = [-4]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, 0.002)

    fillet_indices = [5]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, 0.001)


    # Gate wires 
    for i in range(1, floor(p.CathodeRadialPosition/p.BottomStackWireSpacing)+1):
      # First center coordinates then a point the circle passes
      Sketch.addCircle(p.BottomStackFirstWireRadialPosition+p.BottomStackWireSpacing/4-i*p.BottomStackWireSpacing, 
                       p.CathodeVerticalPosition * p.ShrinkageFactor,
                       p.BottomStackFirstWireRadialPosition+p.BottomStackWireSpacing/4-i*p.BottomStackWireSpacing, 
                       p.CathodeVerticalPosition * p.ShrinkageFactor+p.CathodeWireDiameter/2)

    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        Cathode = model.addFace(part_doc, [Sketch.result()])
        Cathode.setName("Cathode")
        return Cathode
      else:
        return Sketch

def build_topScreen(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("TopScreen")
    else:
      dont_return = True
    x0 = p.TopStackRadialPosition
    y0 = p.TopScreenVerticalPosition - p.TopScreenHeight
    W  = p.TopScreenWidth
    H  = p.TopScreenHeight

    pts = [
        (x0,     y0),       # bottom-left
        (x0,     y0 + H),   # top-left
        (x0 + W, y0 + H),   # top-right
        (x0 + W, y0),       # bottom-right
    ]

    lines = create_polygon_from_corners(Sketch, pts)
    # It wouldnt anchor with 2 for some reason
    anchors = make_anchors(Sketch, lines, pts, line_indices=[0, 1, 2, 3]) 

    # Fillet all external corners: indices 0,1,2,3,4
    fillet_indices = [0,1,2,3]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.TopScreenDeburringRadius)
    
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        TopScreen = model.addFace(part_doc, [Sketch.result()])
        TopScreen.setName("TopScreen")
        return TopScreen
      else:
        return Sketch

def build_bottomScreen(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("BottomScreen")
    else:
      dont_return = True
    x0 = p.BottomScreenRadialPosition
    y0 = p.BottomScreenVerticalPosition*p.ShrinkageFactor
    W  = p.BottomScreenWidth
    H  = p.BottomScreenHeight

    pts = [
        (x0,     y0),       # bottom-left
        (x0,     y0 + H),   # top-left
        (x0 + W, y0 + H),   # top-right
        (x0 + W, y0),       # bottom-right
    ]

    lines = create_polygon_from_corners(Sketch, pts)
    # It wouldnt anchor with 2 for some reason
    anchors = make_anchors(Sketch, lines, pts, line_indices=[0, 1, 2, 3]) 

    # Fillet all external corners: indices 0,1,2,3,4
    fillet_indices = [0,3]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.BottomScreenDeburringRadius)
    fillet_indices = [1,2]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.BottomScreenRoundingRadius)
    
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        BottomScreen = model.addFace(part_doc, [Sketch.result()])
        BottomScreen.setName("BottomScreen")
        return BottomScreen
      else:
        return Sketch
  
def build_copper_ring(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("CopperRing")
    else:
      dont_return = True
    
    x0 = p.CopperRingRadialPosition
    y0 = p.CopperRingVerticalPosition - p.CopperRingHeight
    W  = p.CopperRingWidth
    H  = p.CopperRingHeight

    pts = [
        (x0,     y0),       # bottom-left
        (x0,     y0 + H),   # top-left
        (x0 + W, y0 + H),   # top-right
        (x0 + W, y0),       # bottom-right
    ]

    lines = create_polygon_from_corners(Sketch, pts)
    # It wouldnt anchor with 2 for some reason
    anchors = make_anchors(Sketch, lines, pts, line_indices=[0, 1, 2, 3]) 

    # Fillet all external corners: indices 0,1,2,3,4
    fillet_indices = [0,1,2,3]
    # Build point mapping: corner i is the shared vertex = end of line[i-1]
    n = len(lines)
    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, p.CopperRingDeburringRadius)
    
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        CopperRing = model.addFace(part_doc, [Sketch.result()])
        CopperRing.setName("CopperRing")
        return CopperRing
      else:
        return Sketch

def build_field_cage(part_doc, p, makeface, Sketch=None): 
    
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("FieldCage")
    else:
      dont_return = True

    S  = p.ShrinkageFactor
    P  = p.WireVerticalPitch
    x_wire = p.WireRadialPosition
    Ytop   = p.WireTopVerticalPosition * S
    R      = p.WireDiameter / 2.0

    # One vertical step = half pitch (matches your original loops)

    def make_wire_block(base_y, count, dy_wire, Sketch):
      """
      Create a base circle at y = base_y and make 'count' copies
      using Sketch Linear copy.

      This follows the official TUI example:
        addTranslation([circle.results()[1]], pt1.coordinates(), pt2.coordinates(), count)
      """
      # Base circle; 3-arg form so results()[1] is the edge, just like in the doc
      base_circle = Sketch.addCircle(x_wire, base_y, R)

      # Start and end points for the translation vector
      # Start point is at the circle center (same as GUI "SketchCircle_xxx/Center")
      start_pt = Sketch.addPoint(x_wire, base_y)
      end_pt   = Sketch.addPoint(x_wire, base_y + dy_wire)

      # Linear copy – syntax exactly as in the TUI script
      multi = Sketch.addTranslation(
          [base_circle.results()[1]],
          start_pt.coordinates(),
          end_pt.coordinates(),
          count         # total objects including the original
      )

      return base_circle, multi

    # ---- Block 1: first 5 wires (your first loop) ----
    # centers: Ytop, Ytop - 0.5P S, ..., Ytop - 2.0P S
    base_y1 = Ytop
    dy_wire = -P/2.0 * S
    WireCircle_1, WireTrans_1 = make_wire_block(base_y1, 4, dy_wire, Sketch)

    # ---- Block 2: next 65 wires (your second loop) ----
    # original code gives first center = Ytop - 2.5P S
    base_y2 = base_y1 + 4 * dy_wire
    dy_wire = -P * S
    WireCircle_2, WireTrans_2 = make_wire_block(base_y2, 65, dy_wire, Sketch)

    # ---- Block 3: last 2 wires (your third loop) ----
    base_y3 = base_y2 + 64.5 * dy_wire
    dy_wire = -P/2.0 * S
    WireCircle_3, WireTrans_3 = make_wire_block(base_y3, 2, dy_wire, Sketch)
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        FieldCage = model.addFace(part_doc, [Sketch.result()])
        FieldCage.setName("FieldCage")
        return FieldCage
      else:
        return Sketch



def build_field_cage_guard(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("FieldCageGuard")
    else:
      dont_return = True
    

    S  = p.ShrinkageFactor
    P  = p.WireVerticalPitch
    x_wire = p.WireRadialPosition
    Ytop   = p.WireTopVerticalPosition * S
    R      = p.WireDiameter / 2.0
    

    def make_wire_block(base_y, count, dy_wire, Sketch):
      """
      Create a base circle at y = base_y and make 'count' copies
      using Sketch Linear copy.

      This follows the official TUI example:
        addTranslation([circle.results()[1]], pt1.coordinates(), pt2.coordinates(), count)
      """
      # Base circle; 3-arg form so results()[1] is the edge, just like in the doc
      x0 = p.GuardRadialPosition
      y0 = p.GuardTopVerticalPosition * p.ShrinkageFactor - p.GuardHeight/2 
      W  = p.GuardWidth
      H  = p.GuardHeight
      pts = [
          (x0,     y0),       # bottom-left
          (x0,     y0 + H),   # top-left
          (x0 + W, y0 + H),   # top-right
          (x0 + W, y0),       # bottom-right
      ]
      lines = create_polygon_from_corners(Sketch, pts)
      # It wouldnt anchor with 2 for some reason
      anchors = make_anchors(Sketch, lines, pts, line_indices=[0, 1, 2, 3]) 

      # Selecting the new edges is very annoying  so we collect all objects here and later again 
      feat_list = Sketch.features()
      before_ids = set()
      for i in range(feat_list.size()):
          obj = feat_list.object(i)          # ModelAPI_Object
          before_ids.add(obj)

      # Fillet all external corners: indices 0,1,2,3,4
      fillet_indices = [0,1,2,3]
      # Build point mapping: corner i is the shared vertex = end of line[i-1]
      fillets = []  # store fillet features here
      n = len(lines)
      for i in fillet_indices:
          vertex = lines[(i - 1) % n].endPoint()
          # FIXME : This keeps failing with poorly conditioned vertices on the RHS - for now this worsk but need to find out why this works sometimes
          fillets.append(Sketch.setFilletWithRadius(vertex, p.GuardWidth/2.005))#lines.append()

      feat_list = Sketch.features()
      arcs = []

      for i in range(feat_list.size()):
          obj = feat_list.object(i)          # ModelAPI_Object
          if obj in before_ids:
              continue                       # existed before fillets
          feat = ModelAPI.objectToFeature(obj)   # convert to ModelAPI_Feature
          if not feat:
              continue
          if feat.getKind() == "SketchArc":
              arcs.append(SketchAPI_Arc(feat))

      start_pt = Sketch.addPoint(x0 + W / 2, y0 + H /2)
      end_pt   = Sketch.addPoint(x0 + W / 2, y0 + H /2 + dy_wire)
      
      objects = []
      objects.extend(lines) 
      objects.extend(arcs)

      # Linear copy – syntax exactly as in the TUI script
      multi = Sketch.addTranslation(
          objects,
          start_pt.coordinates(),
          end_pt.coordinates(),
          count         # total objects including the original
      )

      return lines, multi

    base_y1 = Ytop
    dy_wire = -P * S
    res = make_wire_block(base_y1, int(p.GuardNumber) , dy_wire, Sketch)
    
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        FieldCageGuard = model.addFace(part_doc, [Sketch.result()])
        FieldCageGuard.setName("FieldCageGuard")
        return FieldCageGuard
      else:
        return Sketch

def build_bottom_pmts(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("BottomPMT")
    else:
      dont_return = True
    # Not Going to try to keep track we will just take the difference again
    feat_list = Sketch.features()
    before_ids = set()
    for i in range(feat_list.size()):
        obj = feat_list.object(i) 
        before_ids.add(obj)

    # ------------------------------------------------------------------
    # Convenience aliases
    # ------------------------------------------------------------------
    cx = p.BottomPMTsArrayPMTFirstRadialPosition
    QD = p.PMTQuartzDiameter
    QH = p.PMTQuartzHeight
    BD = p.PMTBodyDiameter
    BH = p.PMTBodyHeight
    BaD = p.PMTBasisDiameter
    BaH = p.PMTBasisHeight
    R  = QD / 2.0
    y_ref = (p.BottomPMTsArrayReflectorVerticalPosition
             - p.BottomPMTsArrayReflectorPMTHoleDepth) * p.ShrinkageFactor
    # top of quartz cylinder
    y_q_top = y_ref
    # bottom of quartz cylinder
    y_q_bot = y_q_top - QH
    # circle center (as you gave)
    cy = (p.BottomPMTsArrayReflectorVerticalPosition * p.ShrinkageFactor
          - p.BottomPMTsArrayReflectorPMTHoleDepth * p.ShrinkageFactor
          - p.PMTQuartzHeight)
    # this is equal to y_q_bot, but we keep the expression to mirror COMSOL
    # body and basis levels
    y_body_top    = y_q_bot
    y_body_bottom = y_body_top - BH
    y_basis_bottom = y_body_bottom - BaH
    # x–positions
    x_qL    = cx - QD/2.0           # quartz outer
    x_qR    = cx + QD/2.0
    x_bodyL = cx - BD/2.0           # body outer
    x_bodyR = cx + BD/2.0
    x_basisL = cx - BaD/2.0         # basis outer
    x_basisR = cx + BaD/2.0
    # ------------------------------------------------------------------
    # Intersections of the quartz-circle with body verticals
    #  Circle: (x-cx)^2 + (y-cy)^2 = R^2
    #  Lines:  x = x_bodyL / x_bodyR
    #  We want the LOWER intersection (bell below center): y = cy - dy
    # ------------------------------------------------------------------
    dxL = x_bodyL - cx
    dyL = math.sqrt(max(R*R - dxL*dxL, 0.0))
    y_int_L = cy - dyL

    dxR = x_bodyR - cx
    dyR = math.sqrt(max(R*R - dxR*dxR, 0.0))
    y_int_R = cy - dyR

    left_int_pt  = (x_bodyL, y_int_L)
    right_int_pt = (x_bodyR, y_int_R)

    # Quartz-cylinder / circle contact points (bottom of quartz)
    left_q_bottom  = (x_qL, y_q_bot)
    right_q_bottom = (x_qR, y_q_bot)
    # ------------------------------------------------------------------
    # Build boundary in order (counter-clockwise):
    #
    #  0  left quartz top          (x_qL,   y_q_top)
    #  1  right quartz top         (x_qR,   y_q_top)
    #  2  right quartz bottom      (x_qR,   y_q_bot)
    #  3  right body/circle intercept (x_bodyR, y_int_R)
    #  4  right body bottom        (x_bodyR, y_body_bottom)
    #  5  basis right top          (x_basisR, y_body_bottom)
    #  6  basis right bottom       (x_basisR, y_basis_bottom)
    #  7  basis left bottom        (x_basisL, y_basis_bottom)
    #  8  basis left top           (x_basisL, y_body_bottom)
    #  9  left body bottom         (x_bodyL, y_body_bottom)
    # 10  left body/circle intercept (x_bodyL, y_int_L)
    # 11  left quartz bottom       (x_qL,   y_q_bot)
    # back to 0.
    # ------------------------------------------------------------------
    # 0) Top horizontal (left quartz top -> right quartz top)
    line_top = Sketch.addLine(x_qL, y_q_top, x_qR, y_q_top)
    # 1) Quartz right vertical
    line_q_right = Sketch.addLine(x_qR, y_q_top, x_qR, y_q_bot)
    # 2) Arc from quartz right bottom to body right intercept
    arc_side_R = Sketch.addArc(
        cx, cy,
        right_q_bottom[0], right_q_bottom[1],   # start
        right_int_pt[0],   right_int_pt[1],     # end
        True                                    # flip to True if wrong side
    )
    # 3) Body right vertical (intercept down to body bottom)
    line_body_right = Sketch.addLine(x_bodyR, y_int_R, x_bodyR, y_body_bottom)
    # 4) Basis right top (body right bottom -> basis right outer)
    line_basis_top_R = Sketch.addLine(x_bodyR,  y_body_bottom,
                                      x_basisR, y_body_bottom)
    # 5) Basis right vertical (top -> bottom)
    line_basis_right = Sketch.addLine(x_basisR, y_body_bottom,
                                      x_basisR, y_basis_bottom)
    # 6) Basis bottom (right -> left)
    line_basis_bottom = Sketch.addLine(x_basisR, y_basis_bottom,
                                       x_basisL, y_basis_bottom)
    # 7) Basis left vertical (bottom -> top)
    line_basis_left = Sketch.addLine(x_basisL, y_basis_bottom,
                                     x_basisL, y_body_bottom)
    # 8) Basis left top (basis left top -> body left bottom)
    line_basis_top_L = Sketch.addLine(x_basisL,  y_body_bottom,
                                      x_bodyL,   y_body_bottom)
    # 9) Body left vertical (body bottom -> intercept)
    line_body_left = Sketch.addLine(x_bodyL, y_body_bottom,
                                    x_bodyL, y_int_L)
    # 10) Arc from body left intercept to quartz left bottom
    arc_side_L = Sketch.addArc(
        cx, cy,
        x_bodyL,     y_int_L,    # start
        left_q_bottom[0], left_q_bottom[1],   # end
        True
    )
    # 11) Quartz left vertical (bottom -> top) – closes to top
    line_q_left = Sketch.addLine(x_qL, y_q_bot, x_qL, y_q_top)
    # ------------------------------------------------------------------
    # Coincidence constraints in loop order
    # ------------------------------------------------------------------
    Sketch.setCoincident(line_top.endPoint(),       line_q_right.startPoint())
    Sketch.setCoincident(line_q_right.endPoint(),   arc_side_R.startPoint())
    Sketch.setCoincident(arc_side_R.endPoint(),     line_body_right.startPoint())
    Sketch.setCoincident(line_body_right.endPoint(), line_basis_top_R.startPoint())
    Sketch.setCoincident(line_basis_top_R.endPoint(), line_basis_right.startPoint())
    Sketch.setCoincident(line_basis_right.endPoint(), line_basis_bottom.startPoint())
    Sketch.setCoincident(line_basis_bottom.endPoint(), line_basis_left.startPoint())
    Sketch.setCoincident(line_basis_left.endPoint(), line_basis_top_L.startPoint())
    Sketch.setCoincident(line_basis_top_L.endPoint(), line_body_left.startPoint())
    Sketch.setCoincident(line_body_left.endPoint(),   arc_side_L.startPoint())
    Sketch.setCoincident(arc_side_L.endPoint(),       line_q_left.startPoint())
    Sketch.setCoincident(line_q_left.endPoint(),      line_top.startPoint())

    # --------------   Now we replicate --------------   
    feat_list = Sketch.features()
    arcs = []

    for i in range(feat_list.size()):
        obj = feat_list.object(i)          # ModelAPI_Object
        if obj in before_ids:
            continue                       # existed before fillets
        feat = ModelAPI.objectToFeature(obj)   # convert to ModelAPI_Feature
        if not feat:
            continue
        if (feat.getKind() == "SketchArc"):
            arcs.append(SketchAPI_Arc(feat))
        if  (feat.getKind() == "SketchLine"):
            arcs.append(SketchAPI_Line(feat))

    start_pt = Sketch.addPoint(p.BottomPMTsArrayPMTFirstRadialPosition, (p.BottomPMTsArrayReflectorVerticalPosition-p.BottomPMTsArrayReflectorPMTHoleDepth)*p.ShrinkageFactor)
    end_pt   = Sketch.addPoint(p.BottomPMTsArrayPMTFirstRadialPosition - p.BottomPMTsArrayReflectorPMTHoleRadialPitch, 
                                (p.BottomPMTsArrayReflectorVerticalPosition-p.BottomPMTsArrayReflectorPMTHoleDepth)*p.ShrinkageFactor)

    multi = Sketch.addTranslation(
        arcs,
        start_pt.coordinates(),
        end_pt.coordinates(),
        8         # total objects including the original
    ) 

    # ------------------------------- Need to make half a PMT unfortunately -----------------------------------
    # ----------------- radial positions for axis PMT -------------------
    cx_axis = 0.0                  # CENTER on axis x = 0
    x_q     = QD   / 2.0           # quartz outer radius
    x_body  = BD   / 2.0           # body outer radius
    x_basis = BaD  / 2.0           # basis outer radius

    # Intersection of quartz-circle with body vertical (right side only)
    # Circle: (x - cx_axis)^2 + (y - cy)^2 = R^2,  x = x_body
    dxR = x_body - cx_axis
    dyR = math.sqrt(max(R*R - dxR*dxR, 0.0))
    y_int_R = cy - dyR             # lower intersection (bell below center)

    # ------------------------------------------------------------------
    # Build HALF boundary (x >= 0), counter-clockwise:
    #
    #  A  axis top           (0,       y_q_top)
    #  B  quartz right top   (x_q,     y_q_top)
    #  C  quartz right bottom(x_q,     y_q_bot)
    #  D  body/circle int R  (x_body,  y_int_R)
    #  E  body right bottom  (x_body,  y_body_bottom)
    #  F  basis right top    (x_basis, y_body_bottom)
    #  G  basis right bottom (x_basis, y_basis_bottom)
    #  H  axis bottom        (0,       y_basis_bottom)
    # back to A.
    # ------------------------------------------------------------------

    # Axis points and line
    axis_top_pt    = Sketch.addPoint(0.0, y_q_top)
    axis_bottom_pt = Sketch.addPoint(0.0, y_basis_bottom)
    axis_line = Sketch.addLine(0.0, y_basis_bottom, 0.0, y_q_top)
    axis_line.setAuxiliary(False)   # optional: axis as construction line

    # 0) Top horizontal (axis -> quartz outer)
    line_top = Sketch.addLine(0.0,   y_q_top,
                              x_q,   y_q_top)

    # 1) Quartz right vertical
    line_q_right = Sketch.addLine(x_q, y_q_top,
                                  x_q, y_q_bot)

    # 2) Arc from quartz right bottom to body right intercept
    arc_side_R = Sketch.addArc(
        cx_axis, cy,
        x_q,    y_q_bot,   # start
        x_body, y_int_R,   # end
        True               # same orientation as your original
    )

    # 3) Body right vertical (intercept down to body bottom)
    line_body_right = Sketch.addLine(x_body, y_int_R,
                                     x_body, y_body_bottom)

    # 4) Basis right top (body right bottom -> basis right outer)
    line_basis_top_R = Sketch.addLine(x_body,  y_body_bottom,
                                      x_basis, y_body_bottom)

    # 5) Basis right vertical (top -> bottom)
    line_basis_right = Sketch.addLine(x_basis, y_body_bottom,
                                      x_basis, y_basis_bottom)

    # 6) Basis bottom (right -> axis)
    line_basis_bottom = Sketch.addLine(x_basis,      y_basis_bottom,
                                       0.0,          y_basis_bottom)

    # ------------------------------------------------------------------
    # Coincidence constraints to close the loop
    # ------------------------------------------------------------------
    Sketch.setCoincident(axis_line.endPoint(),     line_top.startPoint())
    Sketch.setCoincident(line_top.endPoint(),      line_q_right.startPoint())
    Sketch.setCoincident(line_q_right.endPoint(),  arc_side_R.startPoint())
    Sketch.setCoincident(arc_side_R.endPoint(),    line_body_right.startPoint())
    Sketch.setCoincident(line_body_right.endPoint(), line_basis_top_R.startPoint())
    Sketch.setCoincident(line_basis_top_R.endPoint(), line_basis_right.startPoint())
    Sketch.setCoincident(line_basis_right.endPoint(), line_basis_bottom.startPoint())
    Sketch.setCoincident(line_basis_bottom.endPoint(), axis_line.startPoint())


    if MODELDO: model.do()

    if not dont_return:
      if makeface:
          PMTFace = model.addFace(part_doc, [Sketch.result()])
          PMTFace.setName("BottomPMT")
          return PMTFace
      else:
          return Sketch

def build_top_pmts(part_doc, p, makeface, Sketch=None):
    """

    TODO This needs major simplification - but not right now Might be best to merge this with above function and simply have one that produces a rotatable PMT that is then positioned

                                    PMTConnectionPoint
                                        ________
                        line_basis_left |      | line_basis_right
                    line_basis_top_L  --       --  line_basis_top_L
                  line_body_left     |           |     line_body_right
                                     |           |     
                                     |           |     
                                     |           |     
                                    _|           |_
                  arc_side_L      /                 \        arc_side_R
      PMTWindowConnectionLeft    |                   |     PMTWindowConnectionRight
                                 |___________________|
                                        PMTWindow
    
    """
    
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("TopPMT")
    else:
      dont_return = True

    # Not going to try to keep track manually – take the difference later to know what is in this structure
    feat_list = Sketch.features()
    before_ids = set()
    for i in range(feat_list.size()):
        obj = feat_list.object(i)
        before_ids.add(obj)

    # ------------------------------------------------------------------
    # Convenience aliases – mirror of bottom, but with Top* parameters
    # ------------------------------------------------------------------
    cx = p.TopPMTsArrayPMTFirstRadialPosition

    QD = p.PMTQuartzDiameter
    QH = p.PMTQuartzHeight
    BD = p.PMTBodyDiameter
    BH = p.PMTBodyHeight
    BaD = p.PMTBasisDiameter
    BaH = p.PMTBasisHeight

    R = QD / 2.0

    # For bottom: y_ref = (ReflectorZ - HoleDepth)*Shrinkage
    # For top (hole above reflector), mirror with +HoleDepth (no shrink unless you have it):
    y_ref = p.TopPMTsArrayReflectorVerticalPosition + p.TopPMTsArrayReflectorPMTHoleDepth

    # bottom of quartz cylinder (inside hole)
    y_q_bot = y_ref
    # top of quartz cylinder
    y_q_top = y_q_bot + QH

    # Circle center for dome pointing UP.
    # Bottom used: cy = ReflectorZ*Shrink - HoleDepth*Shrink - QH = y_q_bot
    # Here we mirror: cy above cylinder bottom by QH:
    cy = y_q_bot + QH  # adjust to your COMSOL expression if needed

    
    # x–positions (same logic as bottom)
    x_qL = cx - QD / 2.0           # quartz outer
    x_qR = cx + QD / 2.0
    x_bodyL = cx - BD / 2.0        # body outer
    x_bodyR = cx + BD / 2.0
    x_basisL = cx - BaD / 2.0      # basis outer
    x_basisR = cx + BaD / 2.0


    # body and basis levels ABOVE the quartz (PMT “upside down” compared to bottom)
    # 1) First compute intercepts from the circle (unchanged)
    dxL = x_bodyL - cx
    dyL = math.sqrt(max(R*R - dxL*dxL, 0.0))
    y_int_L = cy + dyL

    dxR = x_bodyR - cx
    dyR = math.sqrt(max(R*R - dxR*dxR, 0.0))
    y_int_R = cy + dyR

    # 2) Use the intercept (bell–body joint) as the *top* of the body
    #    (for a symmetric PMT you expect y_int_L ≈ y_int_R, but keep both)
    y_body_top_L = y_int_L
    y_body_top_R = y_int_R

    # 3) Body extends *down* from the intercept by BH
    y_body_bottom_L = y_body_top_L + BH
    y_body_bottom_R = y_body_top_R + BH

    # 4) Basis extends further down by BaH
    y_basis_bottom_L = y_body_bottom_L + BaH
    y_basis_bottom_R = y_body_bottom_R + BaH


    # ------------------------------------------------------------------
    # Intersections of the quartz-circle with body verticals
    #  Circle: (x - cx)^2 + (y - cy)^2 = R^2
    #  Lines:  x = x_bodyL / x_bodyR
    #  Now dome is ABOVE center → we want the UPPER intersection: y = cy + dy
    # ------------------------------------------------------------------
    dxL = x_bodyL - cx
    dyL = math.sqrt(max(R * R - dxL * dxL, 0.0))
    y_int_L = cy + dyL

    dxR = x_bodyR - cx
    dyR = math.sqrt(max(R * R - dxR * dxR, 0.0))
    y_int_R = cy + dyR

    left_int_pt = (x_bodyL, y_int_L)
    right_int_pt = (x_bodyR, y_int_R)

    # Quartz-cylinder / circle contact points (TOP of quartz, dome up)
    left_q_top = (x_qL, y_q_top)
    right_q_top = (x_qR, y_q_top)

    # ------------------------------------------------------------------
    # Build boundary in order (counter-clockwise, PMT “upside down”):
    #
    #  0  right quartz bottom       (x_qR,    y_q_bot)
    #  1  left quartz bottom        (x_qL,    y_q_bot)
    #  2  left quartz top           (x_qL,    y_q_top)
    #  3  left body/circle intercept (x_bodyL, y_int_L)
    #  4  left body top             (x_bodyL, y_body_top)
    #  5  basis left top            (x_basisL, y_body_top)
    #  6  basis left bottom         (x_basisL, y_basis_bottom)
    #  7  basis right bottom        (x_basisR, y_basis_bottom)
    #  8  basis right top           (x_basisR, y_body_top)
    #  9  right body top            (x_bodyR, y_body_top)
    # 10  right body/circle intercept (x_bodyR, y_int_R)
    # 11  right quartz top          (x_qR,   y_q_top)
    # back to 0.
    #
    # We keep naming/structure analogous to the bottom function:
    # two arcs joining the quartz-top to the body intercepts,
    # but now the arcs are ABOVE instead of below.
    # ------------------------------------------------------------------

    # Basis bottom (left -> right)
    PMTConnectionPoint = Sketch.addLine(x_basisL, y_basis_bottom_L,
                                      x_basisR, y_basis_bottom_R)

    # Basis right vertical (bottom -> top)
    line_basis_right = Sketch.addLine(x_basisR, y_basis_bottom_R,
                                      x_basisR, y_body_bottom_R)

    # Basis right top (basis right top -> body right bottom)
    line_basis_top_R = Sketch.addLine(x_basisR, y_body_bottom_R,
                                      x_bodyR,  y_body_bottom_R)

    # Body right vertical (body bottom -> intercept)
    line_body_right = Sketch.addLine(x_bodyR, y_body_bottom_R,
                                    x_bodyR, y_int_R)

    # Arc from body right intercept to quartz right top
    arc_side_R = Sketch.addArc(
        cx, cy,
        x_bodyR,    y_int_R,     # start on body
        right_q_top[0], right_q_top[1],  # end on quartz
        True
    )

    # Quartz right vertical (top -> bottom)
    PMTWindowConnectionRight = Sketch.addLine(x_qR, y_q_top, x_qR, y_q_bot)

    # Quartz bottom (right -> left)
    PMTWindow = Sketch.addLine(x_qR, y_q_bot, x_qL, y_q_bot)

    # Quartz left vertical (bottom -> top)
    PMTWindowConnectionLeft = Sketch.addLine(x_qL, y_q_bot, x_qL, y_q_top)

    # Arc from quartz left top to body left intercept
    arc_side_L = Sketch.addArc(
        cx, cy,
        left_q_top[0], left_q_top[1],   # start on quartz
        x_bodyL,      y_int_L,          # end on body
        True
    )

    # Body left vertical (intercept -> body bottom)
    line_body_left = Sketch.addLine(x_bodyL, y_int_L,
                                    x_bodyL, y_body_bottom_L)

    # Basis left top (body left bottom -> basis left top)
    line_basis_top_L = Sketch.addLine(x_bodyL,      y_body_bottom_L,
                                      x_basisL,     y_body_bottom_L)

    # Basis left vertical (top -> bottom)
    line_basis_left = Sketch.addLine(x_basisL, y_body_bottom_L,
                                    x_basisL, y_basis_bottom_L)

    # ------------------------------------------------------------------
    # Coincidence constraints in loop order
    # ------------------------------------------------------------------
    Sketch.setCoincident(PMTWindow.endPoint(),   PMTWindowConnectionLeft.startPoint())
    Sketch.setCoincident(PMTWindowConnectionLeft.endPoint(),     arc_side_L.startPoint())
    Sketch.setCoincident(arc_side_L.endPoint(),      line_body_left.startPoint())
    Sketch.setCoincident(line_body_left.endPoint(),  line_basis_top_L.startPoint())
    Sketch.setCoincident(line_basis_top_L.endPoint(), line_basis_left.startPoint())
    Sketch.setCoincident(line_basis_left.endPoint(), PMTConnectionPoint.startPoint())
    Sketch.setCoincident(PMTConnectionPoint.endPoint(), line_basis_right.startPoint())
    Sketch.setCoincident(line_basis_right.endPoint(),  line_basis_top_R.startPoint())
    Sketch.setCoincident(line_basis_top_R.endPoint(),  line_body_right.startPoint())
    Sketch.setCoincident(line_body_right.endPoint(),   arc_side_R.startPoint())
    Sketch.setCoincident(arc_side_R.endPoint(),        PMTWindowConnectionRight.startPoint())
    Sketch.setCoincident(PMTWindowConnectionRight.endPoint(),      PMTWindow.startPoint())


    # --------------   Now we replicate (exact same pattern) --------------   
    feat_list = Sketch.features()
    arcs_and_lines = []

    for i in range(feat_list.size()):
        obj = feat_list.object(i)          # ModelAPI_Object
        if obj in before_ids:
            continue                       # existed before
        feat = ModelAPI.objectToFeature(obj)        # convert to ModelAPI_Feature
        if not feat:
            continue
        if feat.getKind() == "SketchArc":
            arcs_and_lines.append(SketchAPI_Arc(feat))
        if feat.getKind() == "SketchLine":
            arcs_and_lines.append(SketchAPI_Line(feat))

    # Translation reference: same style as bottom, but with Top* params
    start_pt = Sketch.addPoint(
        p.TopPMTsArrayPMTFirstRadialPosition,
        y_q_bot  # any convenient y on the quartz cylinder; here its bottom
    )

    end_pt = Sketch.addPoint(
        p.TopPMTsArrayPMTFirstRadialPosition - p.TopPMTsArrayReflectorPMTHoleRadialPitch,
        y_q_bot
    )

    multi = Sketch.addTranslation(
        arcs_and_lines,
        start_pt.coordinates(),
        end_pt.coordinates(),
        8    # total objects including the original
    )

    # ---------------------- Make half a PMT -----------------------------------------------------
    # Axis PMT: center on x = 0, use radii for radial positions
    cx_axis = 0.0

    x_q     = QD  / 2.0   # quartz outer radius
    x_body  = BD  / 2.0   # body outer radius
    x_basis = BaD / 2.0   # basis outer radius

    # Intersections of the quartz-circle with body vertical (right side)
    # Circle: (x - cx_axis)^2 + (y - cy)^2 = R^2,  x = x_body
    # Dome is ABOVE center → take upper intersection: y = cy + dy
    dxR = x_body - cx_axis
    dyR = math.sqrt(max(R*R - dxR*dxR, 0.0))
    y_int_R = cy + dyR

    # Body and basis levels below the dome
    y_body_top    = y_int_R
    y_body_bottom = y_body_top + BH
    y_basis_bottom = y_body_bottom + BaH

    # Quartz–circle contact point (top of quartz, dome up)
    right_q_top = (x_q, y_q_top)

    # ------------------------------------------------------------------
    # Build HALF boundary (x >= 0), counter-clockwise:
    #
    #  A  axis bottom         (0,        y_basis_bottom)
    #  B  axis mid            (0,        y_q_bot)
    #  C  quartz right bottom (x_q,      y_q_bot)
    #  D  quartz right top    (x_q,      y_q_top)
    #  E  body/circle int R   (x_body,   y_int_R)
    #  F  body right bottom   (x_body,   y_body_bottom)
    #  G  basis right top     (x_basis,  y_body_bottom)
    #  H  basis right bottom  (x_basis,  y_basis_bottom)
    #  back to A via basis bottom (x-basis -> 0).
    # ------------------------------------------------------------------

    # Axis line from basis bottom up to quartz bottom
    axis_line = Sketch.addLine(0.0, y_basis_bottom, 0.0, y_q_bot)
    axis_line.setAuxiliary(False)  # if you want it as a construction axis

    # Quartz bottom: axis -> quartz right
    line_q_bottom = Sketch.addLine(0.0,   y_q_bot,
                                   x_q,   y_q_bot)

    # Quartz right vertical: bottom -> top
    line_q_right = Sketch.addLine(x_q, y_q_bot,
                                  x_q, y_q_top)

    # Arc from body right intercept to quartz right top (dome)
    arc_side_R = Sketch.addArc(
        cx_axis, cy,
        x_body,      y_int_R,          # start on body
        right_q_top[0], right_q_top[1],# end on quartz top
        True
    )

    # Body right vertical: from intercept/top down to body bottom
    line_body_right = Sketch.addLine(x_body, y_body_top,
                                     x_body, y_body_bottom)

    # Basis top (body bottom -> basis right)
    line_basis_top_R = Sketch.addLine(x_body,   y_body_bottom,
                                      x_basis,  y_body_bottom)

    # Basis right vertical (top -> bottom)
    line_basis_right = Sketch.addLine(x_basis,  y_body_bottom,
                                      x_basis,  y_basis_bottom)

    # Basis bottom (basis right bottom -> axis bottom)
    line_basis_bottom = Sketch.addLine(x_basis,     y_basis_bottom,
                                       0.0,         y_basis_bottom)

    # ------------------------------------------------------------------
    # Coincidence constraints in loop order
    # ------------------------------------------------------------------
    # Axis ↔ quartz bottom
    Sketch.setCoincident(axis_line.endPoint(),      line_q_bottom.startPoint())
    # Quartz bottom ↔ quartz right vertical
    Sketch.setCoincident(line_q_bottom.endPoint(),  line_q_right.startPoint())
    # Quartz right vertical ↔ arc (quartz top)
    Sketch.setCoincident(line_q_right.endPoint(),   arc_side_R.endPoint())
    # Arc ↔ body right vertical (body top)
    Sketch.setCoincident(arc_side_R.startPoint(),   line_body_right.startPoint())
    # Body right vertical ↔ basis top
    Sketch.setCoincident(line_body_right.endPoint(), line_basis_top_R.startPoint())
    # Basis top ↔ basis right vertical
    Sketch.setCoincident(line_basis_top_R.endPoint(), line_basis_right.startPoint())
    # Basis right vertical ↔ basis bottom
    Sketch.setCoincident(line_basis_right.endPoint(), line_basis_bottom.startPoint())
    # Basis bottom ↔ axis bottom (close loop)
    Sketch.setCoincident(line_basis_bottom.endPoint(), axis_line.startPoint())

    if MODELDO: model.do()

    if not dont_return:
      if makeface:
          PMTFace = model.addFace(part_doc, [Sketch.result()])
          PMTFace.setName("TopPMT")
          return PMTFace
      else:
          return Sketch


def build_bell(part_doc, p, makeface, Sketch = None):
    x0 = 0.0
    y_top    = p.BellVerticalPosition
    D_outer  = p.BellOuterDiameter
    H_bell   = p.BellHeight
    t_minor  = p.BellMinorThickness
    t_major  = p.BellMajorThickness
    y_thick  = p.BellThicknessIncreaseVerticalPosition

    pts = [
        (x0,                     y_top),                      # 0 top-left
        (D_outer,                y_top),                      # 1 top-right
        (D_outer,                y_top - H_bell),             # 2 right-bottom thin
        (D_outer - t_minor,      y_top - H_bell),             # 3 inner at thin
        (D_outer - t_minor,      y_thick),                    # 4 inner up to thick transition
        (D_outer - t_major,      y_thick),                    # 5 inner at thick
        (D_outer - t_major,      y_top - t_major),            # 6 inner top of thick section
        (x0,                     y_top - t_major),            # 7 left inner
    ]
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("Bell")
    else:
      dont_return = True
    lines = create_polygon_from_corners(Sketch, pts)
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        Bell = model.addFace(part_doc, [Sketch.result()])
        Bell.setName("Bell")
        return Bell
      else:
        return Sketch



def build_gate_insulating_frame(part_doc, p, makeface, Sketch=None):
    pts = [
        [p.TopStackInsulationRadialPosition, p.GateInsulationVerticalPosition],
        [p.TopStackInsulationRadialPosition, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG+p.GateInsulationDimensionF],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE+p.GateInsulationDimensionD, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG+p.GateInsulationDimensionF],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE+p.GateInsulationDimensionD, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG+p.GateInsulationDimensionF-p.GateInsulationDimensionB],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE+p.GateInsulationDimensionD+p.GateInsulationDimensionC, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG+p.GateInsulationDimensionF-p.GateInsulationDimensionB],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE+p.GateInsulationDimensionD+p.GateInsulationDimensionC, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG+p.GateInsulationDimensionF-p.GateInsulationDimensionB-p.GateInsulationDimensionA],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE+p.GateInsulationDimensionD+p.GateInsulationDimensionC-p.GateInsulationDimensionK, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG+p.GateInsulationDimensionF-p.GateInsulationDimensionB-p.GateInsulationDimensionA],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE+p.GateInsulationDimensionD+p.GateInsulationDimensionC-p.GateInsulationDimensionK, p.GateInsulationVerticalPosition+p.GateInsulationDimensionG+p.GateInsulationDimensionF-p.GateInsulationDimensionB-p.GateInsulationDimensionA+p.GateInsulationDimensionJ],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE, p.GateInsulationVerticalPosition-p.GateInsulationDimensionH],
        [p.TopStackInsulationRadialPosition+p.GateInsulationDimensionE, p.GateInsulationVerticalPosition],
    ]
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("GateInsulatingFrame")
    else:
      dont_return = True
    
    lines = create_polygon_from_corners(Sketch, pts)
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        GateInsulatingFrame = model.addFace(part_doc, [Sketch.result()])
        GateInsulatingFrame.setName("GateInsulatingFrame")
        return GateInsulatingFrame
      else:
        return Sketch

def build_anode_insulating_frame(part_doc, p, makeface, Sketch=None):
    pts = [
        [p.TopStackInsulationRadialPosition, p.AnodeInsulationVerticalPosition],
        [p.TopStackInsulationRadialPosition, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA],
        [p.TopStackInsulationRadialPosition, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionN],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionK, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionN],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionK, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD + p.AnodeInsulationDimensionE, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD + p.AnodeInsulationDimensionE, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF - p.AnodeInsulationDimensionG],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD + p.AnodeInsulationDimensionE - p.AnodeInsulationDimensionH, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF - p.AnodeInsulationDimensionG],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD + p.AnodeInsulationDimensionE - p.AnodeInsulationDimensionH, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF - p.AnodeInsulationDimensionG + p.AnodeInsulationDimensionL],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD + p.AnodeInsulationDimensionE - p.AnodeInsulationDimensionH - p.AnodeInsulationDimensionI, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF - p.AnodeInsulationDimensionG + p.AnodeInsulationDimensionL],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionC + p.AnodeInsulationDimensionD + p.AnodeInsulationDimensionE - p.AnodeInsulationDimensionH - p.AnodeInsulationDimensionI, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF - p.AnodeInsulationDimensionG + p.AnodeInsulationDimensionL + p.AnodeInsulationDimensionJ],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionK, p.AnodeInsulationVerticalPosition + p.AnodeInsulationDimensionA + p.AnodeInsulationDimensionB - p.AnodeInsulationDimensionF - p.AnodeInsulationDimensionG + p.AnodeInsulationDimensionL + p.AnodeInsulationDimensionJ],
        [p.TopStackInsulationRadialPosition + p.AnodeInsulationDimensionK, p.AnodeInsulationVerticalPosition],
    ]
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("AnodeInsulatingFrame")
    else:
      dont_return = True
    
    lines = create_polygon_from_corners(Sketch, pts)
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        AnodeInsulatingFrame = model.addFace(part_doc, [Sketch.result()])
        AnodeInsulatingFrame.setName("AnodeInsulatingFrame")
        return AnodeInsulatingFrame
      else:
        return Sketch

def build_cathode_insulating_frame(part_doc, p, makeface, Sketch=None):
    pts = [
        [p.WireRadialPosition + 0.00114,                                     (p.CathodeInsulationVerticalPosition + 0.0265)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003,                             (p.CathodeInsulationVerticalPosition + 0.0265)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003,                             (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275,                    (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275,                    (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275 - 0.002,            (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275 - 0.002,            (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013 + 0.0055)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275 - 0.002 - 0.002,    (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013 + 0.0055)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275 - 0.002 - 0.002,    (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013 + 0.0055 + 0.0045)*p.ShrinkageFactor + 0.00036],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275 - 0.002 - 0.002 - 0.018, (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013 + 0.0055 + 0.0045)*p.ShrinkageFactor + 0.00036],
        [p.WireRadialPosition + 0.00114 + 0.003 + 0.0275 - 0.002 - 0.002 - 0.018, (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013 + 0.0055)*p.ShrinkageFactor],
        [p.WireRadialPosition + 0.00114,                                     (p.CathodeInsulationVerticalPosition + 0.0265 - 0.019 - 0.013 + 0.0055)*p.ShrinkageFactor]
    ]
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("CathodeInsulatingFrame")
    else:
      dont_return = True

    
    lines = create_polygon_from_corners(Sketch, pts)
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        CathodeInsulatingFrame = model.addFace(part_doc, [Sketch.result()])
        CathodeInsulatingFrame.setName("CathodeInsulatingFrame")
        return CathodeInsulatingFrame
      else:
        return Sketch

def build_topscreen_insulating_frame(part_doc, p, makeface, Sketch = None):
    pts = [
        [p.TopStackInsulationRadialPosition, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF + p.TopScreenInsulationDimensionD, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF + p.TopScreenInsulationDimensionD, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE - p.TopScreenInsulationDimensionC - p.TopScreenInsulationDimensionB - p.TopScreenInsulationDimensionA],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF + p.TopScreenInsulationDimensionD - p.TopScreenInsulationDimensionJ, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE - p.TopScreenInsulationDimensionC - p.TopScreenInsulationDimensionB - p.TopScreenInsulationDimensionA],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF + p.TopScreenInsulationDimensionD - p.TopScreenInsulationDimensionJ, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE - p.TopScreenInsulationDimensionC - p.TopScreenInsulationDimensionB],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF + p.TopScreenInsulationDimensionD - p.TopScreenInsulationDimensionJ - p.TopScreenInsulationDimensionI, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE - p.TopScreenInsulationDimensionC - p.TopScreenInsulationDimensionB],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF + p.TopScreenInsulationDimensionD - p.TopScreenInsulationDimensionJ - p.TopScreenInsulationDimensionI, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE - p.TopScreenInsulationDimensionC],
        [p.TopStackInsulationRadialPosition + p.TopScreenInsulationDimensionF, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE - p.TopScreenInsulationDimensionC],
        [p.TopStackInsulationRadialPosition, p.TopScreenInsulationVerticalPosition + p.TopScreenInsulationDimensionG - p.TopScreenInsulationDimensionE - p.TopScreenInsulationDimensionC]
    ]
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("TopScreenInsulatingFrame")
    else:
      dont_return = True
    lines = create_polygon_from_corners(Sketch, pts)
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        TopScreenInsulatingFrame = model.addFace(part_doc, [Sketch.result()])
        TopScreenInsulatingFrame.setName("TopScreenInsulatingFrame")
        return TopScreenInsulatingFrame
      else:
        return Sketch

def build_bottomstack_insulating_frame(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("BottomStackInsulatingFrame")
    else:
      dont_return = True
    
    """
    we buidl a rectangle
     __
    |  |
    |  |
    
    And then remove circles out of both sides to get the grooves
    Here this is done with edges only

    We are going anticlockwise with line creation
    """
    # Rectangle parameters
    W  = p.BottomStackInsulationWidth
    H  = p.BottomStackInsulationHeight * p.ShrinkageFactor

    x0 = p.BottomStackInsulationRadialPosition
    y0 = (p.BottomStackInsulationVerticalPosition * p.ShrinkageFactor
          - p.BottomStackInsulationHeight * p.ShrinkageFactor)  # bottom of rectangle

    # Top of rectangle (for clarity)
    y_top = y0 + H

    # Circle/groove parameters
    r      = p.BottomStackInsulationGrooveDiameter / 2.0
    pitch  = p.BottomStackInsulationGrooveVerticalPitch * p.ShrinkageFactor

    # Center of the topmost groove row (for both left and right stacks)
    yc_top = (p.BottomStackInsulationVerticalPosition * p.ShrinkageFactor
              - p.BottomStackInsulationGrooveTopVerticalPosition * p.ShrinkageFactor)

    # Number of groove rows
    n_rows = int(p.BottomStackInsulationHeight / p.BottomStackInsulationGrooveVerticalPitch)

    # X-positions of the two groove columns
    x_left  = x0
    x_right = x0 + W

    ### Compute First and final circle centers and top bottom coordinates
    yc_first = yc_top
    yc_final = yc_top - (n_rows - 1) * pitch

    yc_first_top_intersect = yc_first + r
    yc_first_bottom_intersect = yc_first - r

    yc_final_top_intersect = yc_final + r
    yc_final_bottom_intersect = yc_final - r

    # Create the starting arcs for each side (top for left bottom for right since anticlockwise)
    start_left = Sketch.addPoint(x_left, yc_first)
    end_left = Sketch.addPoint(x_left, yc_final)
    first_arc_L = Sketch.addArc(
        x_left, yc_first,
        x_left, yc_first_top_intersect,    # start point
        x_left, yc_first_bottom_intersect, # end point
        True
    )
    start_right = Sketch.addPoint(x_right, yc_first)
    end_right = Sketch.addPoint(x_right, yc_final)
    last_arc_R = Sketch.addArc(
        x_right, yc_final,
        x_right, yc_final_bottom_intersect,    # start point
        x_right, yc_final_top_intersect, # end point
        True
    )

    left_translation = Sketch.addTranslation(
        [first_arc_L],
        start_left.coordinates(),
        end_left.coordinates(),
        n_rows,         # total objects including the original
        True       # Place within range not repeat range
    )
    arcs_left = [first_arc_L] + [SketchAPI_Arc(i) for i in left_translation.translated()]
    right_translation = Sketch.addTranslation(
        [last_arc_R],
        end_right.coordinates(),
        start_right.coordinates(),
        n_rows,         # total objects including the original
        True       # Place within range not repeat range
    )
    arcs_right = [last_arc_R] + [SketchAPI_Arc(i) for i in right_translation.translated()]

    # Now we have a list of arcs and we connect them anticlockwise for each column first - so we get all the internal lines (not the top rectangle and bottom rectangle ones yet)
    lines_left = []
    for i in range(1, len(arcs_left)):
      lines_left.append(
        Sketch.addLine(arcs_left[i-1].endPoint().x(), arcs_left[i-1].endPoint().y(),
                       arcs_left[i].startPoint().x(), arcs_left[i].startPoint().y())
                       )
      # And make everything coincident
      Sketch.setCoincident(arcs_left[i-1].endPoint(), lines_left[-1].startPoint())
      Sketch.setCoincident(arcs_left[i].startPoint(), lines_left[-1].endPoint())
      
    lines_right = []
    for i in range(1, len(arcs_right)):
      lines_right.append(
        Sketch.addLine(arcs_right[i-1].endPoint().x(), arcs_right[i-1].endPoint().y(),
                       arcs_right[i].startPoint().x(), arcs_right[i].startPoint().y())
                       )
      # And make everything coincident
      Sketch.setCoincident(arcs_right[i-1].endPoint(), lines_right[-1].startPoint())
      Sketch.setCoincident(arcs_right[i].startPoint(), lines_right[-1].endPoint())
      
    # Now we just need to finish up the rectangles, starting with the bottom connecting arcs_left[-1] to corner to corner to arcs_right[0]
    # We append to lines_misc as ordering this is too annoying
    lines_misc = []
    lines_misc.append(
      Sketch.addLine(arcs_left[-1].endPoint().x(), arcs_left[-1].endPoint().y(),
                     p.BottomStackInsulationRadialPosition, (p.BottomStackInsulationVerticalPosition-p.BottomStackInsulationHeight)*p.ShrinkageFactor  
                    )
    )
    Sketch.setCoincident(arcs_left[-1].endPoint(), lines_misc[-1].startPoint())
    lines_misc.append(
      Sketch.addLine(lines_misc[-1].endPoint().x(), lines_misc[-1].endPoint().y(),
                     p.BottomStackInsulationRadialPosition + p.BottomStackInsulationWidth, (p.BottomStackInsulationVerticalPosition-p.BottomStackInsulationHeight)*p.ShrinkageFactor  
                    )
    )
    Sketch.setCoincident(lines_misc[-2].endPoint(), lines_misc[-1].startPoint())
    # Last line connecting to arcs_right
    lines_misc.append(
      Sketch.addLine(lines_misc[-1].endPoint().x(), lines_misc[-1].endPoint().y(),
                     arcs_right[0].startPoint().x(), arcs_right[0].startPoint().y()
                    )
    )
    Sketch.setCoincident(lines_misc[-2].endPoint(), lines_misc[-1].startPoint())
    # And make coincident with the first arc
    Sketch.setCoincident(lines_misc[-1].endPoint(), arcs_right[0].startPoint())

    # Now the top part same logic start at arcs_right -1 go edge egde arcs_left 0 
    lines_misc.append(
      Sketch.addLine(arcs_right[-1].endPoint().x(), arcs_right[-1].endPoint().y(),
                     p.BottomStackInsulationRadialPosition + p.BottomStackInsulationWidth, (p.BottomStackInsulationVerticalPosition)*p.ShrinkageFactor  
                    )
    )
    Sketch.setCoincident(arcs_right[-1].endPoint(), lines_misc[-1].startPoint())
    lines_misc.append(
      Sketch.addLine(lines_misc[-1].endPoint().x(), lines_misc[-1].endPoint().y(),
                     p.BottomStackInsulationRadialPosition, (p.BottomStackInsulationVerticalPosition)*p.ShrinkageFactor  
                    )
    )
    Sketch.setCoincident(lines_misc[-2].endPoint(), lines_misc[-1].startPoint())
    # Last line connecting to arcs_right
    lines_misc.append(
      Sketch.addLine(lines_misc[-1].endPoint().x(), lines_misc[-1].endPoint().y(),
                     arcs_left[0].startPoint().x(), arcs_left[0].startPoint().y()
                    )
    )
    Sketch.setCoincident(lines_misc[-2].endPoint(), lines_misc[-1].startPoint())
    # And make coincident with the first arc
    Sketch.setCoincident(lines_misc[-1].endPoint(), arcs_left[0].startPoint())

    if MODELDO: model.do()
    if makeface:
      BottomStackInsulatingFrame = model.addFace(part_doc, [Sketch.result()])
      BottomStackInsulatingFrame.setName("BottomStackInsulatingFrame")
      return BottomStackInsulatingFrame
    else:
      return Sketch

def build_bottom_PMT_reflectors(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("BottomPMTReflectors")
    else:
      dont_return = True

    # In comsol this was defined as left offset to provide a pin the bottom screening sits on, that part will be made last so we can abuse the repeating architecture
    dx = -p.BottomPMTsArrayReflectorPMTHoleRadialPitch
    poly_left = [
        [
            p.BottomPMTsArrayPMTFirstRadialPosition - p.BottomPMTsArrayReflectorPMTHoleHousingDiameter/2,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorHeight*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition - p.BottomPMTsArrayReflectorPMTHoleHousingDiameter/2,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorPMTHoleDepth*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition - p.BottomPMTsArrayReflectorPMTHoleDiameter/2 + p.BottomPMTsArrayReflectorPMTHoleChamfer,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorPMTHoleDepth*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition - p.BottomPMTsArrayReflectorPMTHoleDiameter/2 + p.BottomPMTsArrayReflectorPMTHoleChamfer,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorPMTHoleChamfer*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition - p.BottomPMTsArrayReflectorPMTHoleDiameter/2,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
        ],
    ]

    poly_right = [
        [
            p.BottomPMTsArrayPMTFirstRadialPosition + p.BottomPMTsArrayReflectorPMTHoleDiameter/2 + dx,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition + p.BottomPMTsArrayReflectorPMTHoleDiameter/2 - p.BottomPMTsArrayReflectorPMTHoleChamfer + dx,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorPMTHoleChamfer*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition + p.BottomPMTsArrayReflectorPMTHoleDiameter/2 - p.BottomPMTsArrayReflectorPMTHoleChamfer + dx,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorPMTHoleDepth*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition + p.BottomPMTsArrayReflectorPMTHoleHousingDiameter/2 + dx,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorPMTHoleDepth*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayPMTFirstRadialPosition + p.BottomPMTsArrayReflectorPMTHoleHousingDiameter/2 + dx,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor - p.BottomPMTsArrayReflectorHeight*p.ShrinkageFactor
        ],
    ]

    lines = create_polygon_from_corners(Sketch, poly_left+poly_right)
    anchors = make_anchors(Sketch, lines, poly_left+poly_right, line_indices=[0, 1, 2, 3]) 
    # Now we repeat this 1 times less than pmts
    start = Sketch.addPoint(poly_left[0][0], poly_left[0][1])
    end =   Sketch.addPoint(poly_left[0][0] + dx, poly_left[0][1])
    left_translation = Sketch.addTranslation(
        lines,
        start.coordinates(),
        end.coordinates(),
        8,         # total objects including the original
        False       # Distance for one element
    )

    # And now we make the right bit the bottom shielding sits on 
    # First we remove the shift we made 
    for i in range(len(poly_right)):
      poly_right[i][0] -= dx

    extra_points = [
        [
            p.BottomPMTsArrayReflectorDiameter/2,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
                - p.BottomPMTsArrayReflectorHeight*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayReflectorDiameter/2,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayReflectorNotchRadialPosition
                + p.BottomPMTsArrayReflectorNotchWidth,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayReflectorNotchRadialPosition
                + p.BottomPMTsArrayReflectorNotchWidth,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
                - p.BottomPMTsArrayReflectorNotchDepth*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayReflectorNotchRadialPosition,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
                - p.BottomPMTsArrayReflectorNotchDepth*p.ShrinkageFactor
        ],
        [
            p.BottomPMTsArrayReflectorNotchRadialPosition,
            p.BottomPMTsArrayReflectorVerticalPosition*p.ShrinkageFactor
        ]
    ]
    lines2 = create_polygon_from_corners(Sketch, poly_right+extra_points)
    anchors = make_anchors(Sketch, lines2, poly_right+extra_points, line_indices=[0, 1, 2, 3]) 

    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        BottomPMTReflectors = model.addFace(part_doc, [Sketch.result()])
        BottomPMTReflectors.setName("BottomPMTReflectors")
        return BottomPMTReflectors
      else:
        return Sketch

def build_top_PMT_reflectors(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("TopPMTReflectors")
    else:
      dont_return = True

    # In comsol this was defined as left offset to provide a pin the bottom screening sits on, that part will be made last so we can abuse the repeating architecture
    dx = -p.TopPMTsArrayReflectorPMTHoleRadialPitch
    # Left side of the top reflector cell (no offset)
    poly_left = [
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                - p.TopPMTsArrayReflectorPMTHoleHousingDiameter/2,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorHeight
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                - p.TopPMTsArrayReflectorPMTHoleHousingDiameter/2,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorPMTHoleDepth
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                - p.TopPMTsArrayReflectorPMTHoleDiameter/2
                + p.TopPMTsArrayReflectorPMTHoleChamfer,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorPMTHoleDepth
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                - p.TopPMTsArrayReflectorPMTHoleDiameter/2
                + p.TopPMTsArrayReflectorPMTHoleChamfer,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorPMTHoleChamfer
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                - p.TopPMTsArrayReflectorPMTHoleDiameter/2,
            p.TopPMTsArrayReflectorVerticalPosition
        ]
    ]

    # Right side of the top reflector cell, shifted in x by dx
    poly_right = [
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                + p.TopPMTsArrayReflectorPMTHoleDiameter/2
                + dx,
            p.TopPMTsArrayReflectorVerticalPosition
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                + p.TopPMTsArrayReflectorPMTHoleDiameter/2
                - p.TopPMTsArrayReflectorPMTHoleChamfer
                + dx,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorPMTHoleChamfer
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                + p.TopPMTsArrayReflectorPMTHoleDiameter/2
                - p.TopPMTsArrayReflectorPMTHoleChamfer
                + dx,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorPMTHoleDepth
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                + p.TopPMTsArrayReflectorPMTHoleHousingDiameter/2
                + dx,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorPMTHoleDepth
        ],
        [
            p.TopPMTsArrayPMTFirstRadialPosition
                + p.TopPMTsArrayReflectorPMTHoleHousingDiameter/2
                + dx,
            p.TopPMTsArrayReflectorVerticalPosition
                + p.TopPMTsArrayReflectorHeight
        ]
    ]
    lines = create_polygon_from_corners(Sketch, poly_left+poly_right)
    anchors = make_anchors(Sketch, lines, poly_left+poly_right, line_indices=[0, 1, 2, 3]) 

    # Now we repeat this 1 times less than pmts
    start = Sketch.addPoint(poly_left[0][0], poly_left[0][1])
    end =   Sketch.addPoint(poly_left[0][0] + dx, poly_left[0][1])
    left_translation = Sketch.addTranslation(
        lines,
        start.coordinates(),
        end.coordinates(),
        8,         # total objects including the original
        False      # Distance for one element
    )

    # And now we make the right bit the bottom shielding sits on 
    # First we remove the shift we made 
    for i in range(len(poly_right)):
        poly_right[i][0] -= dx

    extra_points = [
        [p.TopPMTsArrayReflectorDiameter/2, p.TopPMTsArrayReflectorVerticalPosition + p.TopPMTsArrayReflectorHeight],
        [p.TopPMTsArrayReflectorDiameter/2, p.TopPMTsArrayReflectorVerticalPosition],
    ]
    lines2 = create_polygon_from_corners(Sketch, poly_right+extra_points)
    anchors = make_anchors(Sketch, lines2, poly_right+extra_points, line_indices=[0, 1, 2, 3]) 

    if MODELDO: model.do()
    if not dont_return:
      if makeface:
          TopPMTReflectors = model.addFace(part_doc, [Sketch.result()])
          TopPMTReflectors.setName("TopPMTReflectors")
          return TopPMTReflectors
      else:
          return Sketch


def BuildPTFEWall(part_doc, p, makeface, Sketch=None):
    dont_return = False
    if Sketch is None:
      Sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
      Sketch.setName("wallPTFE")
    else:
      dont_return = True
    
    x0 = p.PanelRadialPosition
    y0 = p.PanelVerticalPosition - p.PanelHeight*p.ShrinkageFactor
    W  = p.PanelWidth
    H  = p.PanelHeight*p.ShrinkageFactor

    pts = [
        (x0,     y0),       # bottom-left
        (x0,     y0 + H),   # top-left
        (x0 + W, y0 + H),   # top-right
        (x0 + W, y0),       # bottom-right
    ]

    lines = create_polygon_from_corners(Sketch, pts)
    anchors = make_anchors(Sketch, lines, pts, line_indices=[0, 1, 2, 3]) 
    if MODELDO: model.do()
    if not dont_return:
      if makeface:
        wallPTFE = model.addFace(part_doc, [Sketch.result()])
        wallPTFE.setName("wallPTFE")
        return wallPTFE
      else:
        return Sketch
def face_feature_to_selections(src):
    """
    Turn any 'face-like' object into a list of FACE selections.

    src can be:
      - ModelHighAPI_Selection       (already a selection)
      - feature with .results()      (Cut, Group, etc., possibly multi-volume)
      - feature with .result()       (Face with single result)
      - ModelAPI_Result              (possibly multi-face via subResult)
    """

    # 1) Already a selection?
    if hasattr(src, "type") and callable(getattr(src, "type", None)):
        # Assume it's already a FACE selection; if you need, check src.type() == "FACE"
        return [src]

    # 2) Feature with multiple results (e.g. Cut)
    if hasattr(src, "results"):
        res_list = list(src.results())

    # 3) Feature with single result (e.g. addFace)
    elif hasattr(src, "result"):
        res_list = [src.result()]

    else:
        # 4) Assume it's already a Result
        res_list = [src]

    sels = []
    for r in res_list:
        # r is a Result; check for sub-faces
        if hasattr(r, "numberOfSubs") and r.numberOfSubs() > 0:
            for i in range(r.numberOfSubs()):
                sub_sel = r.subResult(i)  # This is already a ModelHighAPI_Selection
                sels.append(sub_sel)
        else:
            # Single face under this result
            sels.append(model.selection("FACE", r.name()))
    return sels


from collections import defaultdict, deque

def shape_centroid(shape):
    """Centroid of all vertices of a shape."""
    verts = geompy.ExtractShapes(shape, geompy.ShapeType["VERTEX"], True)
    if not verts:
        return None
    sx = sy = sz = 0.0
    for v in verts:
        x, y, z = geompy.PointCoordinates(v)
        sx += x; sy += y; sz += z
    n = len(verts)
    return (sx / n, sy / n, sz / n)

def edge_vertices(edge):
    verts = geompy.ExtractShapes(edge, geompy.ShapeType["VERTEX"], True)
    if len(verts) < 2:
        return None, None
    return verts[0], verts[1]

def analyze_bc_loops(edge_list, bc_name, debug = False):
    """
    Split edge_list into connected components via shared vertices and
    check each component ("loop") separately, using OCC topological
    identity (v.IsSame(u)) to merge coincident vertices.
    """
    if not edge_list:
        print(f"BC '{bc_name}': no edges to analyze.")
        return

    # ------------------------------------------------------------
    # 1) Canonical vertex mapping using IsSame
    # ------------------------------------------------------------
    canonical_vertices = []  # list of unique OCC-vertices
    def canon(v):
        """Return canonical representative for vertex v using IsSame."""
        for u in canonical_vertices:
            if v.IsSame(u):
                return u
        canonical_vertices.append(v)
        return v

    edge_vertices_map = {}           # edge -> (canon_v1, canon_v2)
    vertex_neighbors = defaultdict(set)

    for e in edge_list:
        verts = geompy.ExtractShapes(e, geompy.ShapeType["VERTEX"], True)
        if len(verts) < 2:
            print(f"BC '{bc_name}': edge with <2 vertices; skipping in loop analysis.")
            continue

        v1_raw, v2_raw = verts[0], verts[1]
        v1 = canon(v1_raw)
        v2 = canon(v2_raw)

        edge_vertices_map[e] = (v1, v2)
        vertex_neighbors[v1].add(v2)
        vertex_neighbors[v2].add(v1)

    if not edge_vertices_map:
        print(f"BC '{bc_name}': no usable edges for loop analysis.")
        return

    # ------------------------------------------------------------
    # 2) Connected components in the canonical-vertex graph
    # ------------------------------------------------------------
    unvisited = set(vertex_neighbors.keys())
    components = []

    while unvisited:
        start = next(iter(unvisited))
        queue = deque([start])
        comp_vertices = {start}
        unvisited.remove(start)

        while queue:
            v = queue.popleft()
            for nbr in vertex_neighbors[v]:
                if nbr in unvisited:
                    unvisited.remove(nbr)
                    comp_vertices.add(nbr)
                    queue.append(nbr)
        components.append(comp_vertices)

    # ------------------------------------------------------------
    # 3) Analyze each component => closed/open
    # ------------------------------------------------------------
    for i, comp_vertices in enumerate(components):
        comp_edges = []
        vertex_degree = defaultdict(int)

        for e, (v1, v2) in edge_vertices_map.items():
            if v1 in comp_vertices or v2 in comp_vertices:
                comp_edges.append(e)
                if v1 in comp_vertices:
                    vertex_degree[v1] += 1
                if v2 in comp_vertices:
                    vertex_degree[v2] += 1

        if not comp_edges:
            print(f"BC '{bc_name}', loop {i}: no edges in this component (unexpected).")
            continue

        # Degree-1 canonical vertices are open ends
        open_vertices = [v for v, deg in vertex_degree.items() if deg == 1]

        if open_vertices:
            print(
                f"BC '{bc_name}', loop {i}: OPEN "
                f"({len(open_vertices)} vertex(ces) with degree 1)."
            )
            for j, v in enumerate(open_vertices):
                x, y, z = geompy.PointCoordinates(v)
                print(f"    open vertex {j}: ({x:.6g}, {y:.6g}, {z:.6g})")
        else:
            if debug:
                print(
                    f"BC '{bc_name}', loop {i}: CLOSED "
                    f"({len(comp_edges)} edges, {len(comp_vertices)} vertices)."
                )
def bc_components_from_edges(bc_edges, bc_name):
    """
    Given all edges of a BC shape, split into connected components via vertices.
    Returns a list of components; each component is a set of vertex objects.
    """
    edge_vertices_map = {}
    vertex_neighbors = defaultdict(set)

    for e in bc_edges:
        v1, v2 = edge_vertices(e)
        if v1 is None:
            print(f"BC '{bc_name}': degenerate BC edge; skipping in component build.")
            continue
        edge_vertices_map[e] = (v1, v2)
        vertex_neighbors[v1].add(v2)
        vertex_neighbors[v2].add(v1)

    components = []
    unvisited = set(vertex_neighbors.keys())

    while unvisited:
        start = next(iter(unvisited))
        queue = deque([start])
        comp_vertices = {start}
        unvisited.remove(start)

        while queue:
            v = queue.popleft()
            for nbr in vertex_neighbors[v]:
                if nbr in unvisited:
                    unvisited.remove(nbr)
                    comp_vertices.add(nbr)
                    queue.append(nbr)
        components.append(comp_vertices)

    return components

def component_centroid(comp_vertices):
    """Centroid of a set of vertex objects."""
    if not comp_vertices:
        return None
    sx = sy = sz = 0.0
    for v in comp_vertices:
        x, y, z = geompy.PointCoordinates(v)
        sx += x; sy += y; sz += z
    n = len(comp_vertices)
    return (sx / n, sy / n, sz / n)
def match_component_to_wire(comp_vertices,
                            wire_cache,
                            dist_tol=1e-3,
                            min_frac=0.3):
    """
    Match a BC component (set of vertices) to a *single* domain wire,
    using vertex-based proximity.

    - dist_tol: maximum distance for a BC vertex to count as 'matching'
      some vertex on the wire.
    - min_frac: minimum fraction of BC vertices that must be matched
      for a wire to be accepted.

    Returns: (best_wire, best_index, best_matches, best_avg_dist)
    or (None, None, 0, None) if nothing matches.
    """
    if not comp_vertices:
        return None, None, 0, None

    # Precompute BC vertex coordinates
    bc_coords = []
    for v in comp_vertices:
        bc_coords.append(geompy.PointCoordinates(v))

    best_wire   = None
    best_index  = None
    best_score  = 0
    best_avg    = None

    for w, verts, idx in wire_cache:
        matches = 0
        dsum    = 0.0

        # For each BC vertex, find nearest vertex on this wire
        for (cx, cy, cz) in bc_coords:
            mind = None
            for v_wire, (wx, wy, wz) in verts:
                dx = cx - wx
                dy = cy - wy
                dz = cz - wz
                d = math.sqrt(dx*dx + dy*dy + dz*dz)
                if mind is None or d < mind:
                    mind = d
            if mind is not None and mind < dist_tol:
                matches += 1
                dsum    += mind

        if matches == 0:
            continue

        frac = matches / float(len(bc_coords))
        if frac < min_frac:
            continue

        avg = dsum / float(matches)

        if (matches > best_score) or (matches == best_score and
                                      (best_avg is None or avg < best_avg)):
            best_score = matches
            best_avg   = avg
            best_wire  = w
            best_index = idx

    return best_wire, best_index, best_score, best_avg

def canon_vertex(v, canonical_list):
    """Return canonical representative for vertex v using IsSame."""
    for u in canonical_list:
        if v.IsSame(u):
            return u
    canonical_list.append(v)
    return v

def classify_wire(wire):
    """
    For a given WIRE, classify it based on:
      - open/closed (vertex degree),
      - whether any 'open' vertex lies on x = 0 (within axis_tol).

    Returns: ("closed" | "axis-open" | "other-open"),
            canonical_vertices,
            vertex_neighbors
    """
    verts_raw = geompy.ExtractShapes(wire, geompy.ShapeType["VERTEX"], True)
    if not verts_raw:
        return "other-open", [], {}

    canonical = []
    vertex_neighbors = defaultdict(set)

    # Build edges and neighbor graph on canonical vertices
    edges = geompy.ExtractShapes(wire, geompy.ShapeType["EDGE"], True)
    for e in edges:
        v_raw = geompy.ExtractShapes(e, geompy.ShapeType["VERTEX"], True)
        if len(v_raw) < 2:
            continue
        v1 = canon_vertex(v_raw[0], canonical)
        v2 = canon_vertex(v_raw[1], canonical)
        vertex_neighbors[v1].add(v2)
        vertex_neighbors[v2].add(v1)

    if not vertex_neighbors:
        return "other-open", canonical, vertex_neighbors

    # Degree per canonical vertex
    vertex_degree = {v: len(neigh) for v, neigh in vertex_neighbors.items()}
    open_vertices = [v for v, deg in vertex_degree.items() if deg == 1]

    if not open_vertices:
        loop_type = "closed"
    else:
        # Are all open vertices on the axis?
        all_on_axis = True
        for v in open_vertices:
            x, y, z = geompy.PointCoordinates(v)
            if abs(x) > axis_tol:
                all_on_axis = False
                break
        loop_type = "axis-open" if all_on_axis else "other-open"

    return loop_type, canonical, vertex_neighbors

def build_wire_cache_with_classification(domain_shape):
    """
    Extract all WIREs from domain_shape and classify them.
    Returns a list of dicts:
        {
          "wire": wire_obj,
          "index": idx,
          "type": "closed"|"axis-open"|"other-open",
          "vertices": [(vertex_obj, (x,y,z)), ...]
        }
    """
    wires = geompy.ExtractShapes(domain_shape, geompy.ShapeType["WIRE"], True)
    cache = []

    for idx, w in enumerate(wires):
        loop_type, canon_verts, _ = classify_wire(w)
        vlist = []
        for v in canon_verts:
            x, y, z = geompy.PointCoordinates(v)
            vlist.append((v, (x, y, z)))
        cache.append({
            "wire": w,
            "index": idx,
            "type": loop_type,
            "vertices": vlist,
        })
    return cache

def bc_components_from_edges_with_edges(bc_edges, bc_name):
    """
    Split bc_edges into connected components via shared vertices.
    Returns a list of (component_vertices_set, component_edges_list).
    """
    edge_vertices_map = {}
    vertex_neighbors = defaultdict(set)

    for e in bc_edges:
        verts = geompy.ExtractShapes(e, geompy.ShapeType["VERTEX"], True)
        if len(verts) < 2:
            print(f"BC '{bc_name}': degenerate BC edge; skipping in component build.")
            continue
        v1, v2 = verts[0], verts[1]
        edge_vertices_map[e] = (v1, v2)
        vertex_neighbors[v1].add(v2)
        vertex_neighbors[v2].add(v1)

    components = []
    unvisited = set(vertex_neighbors.keys())

    while unvisited:
        start = next(iter(unvisited))
        queue = deque([start])
        comp_vertices = {start}
        unvisited.remove(start)

        while queue:
            v = queue.popleft()
            for nbr in vertex_neighbors[v]:
                if nbr in unvisited:
                    unvisited.remove(nbr)
                    comp_vertices.add(nbr)
                    queue.append(nbr)

        # Collect edges belonging to this component
        comp_edges = []
        for e, (v1, v2) in edge_vertices_map.items():
            if v1 in comp_vertices or v2 in comp_vertices:
                comp_edges.append(e)

        components.append((comp_vertices, comp_edges))

    return components

def classify_bc_component(comp_vertices, comp_edges):
    """
    Classify a BC component (vertices + edges) as in classify_wire:
      - "closed": no open vertices
      - "axis-open": open vertices exist, all on x=0 (|x|<axis_tol)
      - "other-open": open vertices not all on axis
    """
    if not comp_edges:
        return "other-open"

    canonical = []
    vertex_neighbors = defaultdict(set)

    for e in comp_edges:
        verts = geompy.ExtractShapes(e, geompy.ShapeType["VERTEX"], True)
        if len(verts) < 2:
            continue
        v1 = canon_vertex(verts[0], canonical)
        v2 = canon_vertex(verts[1], canonical)
        vertex_neighbors[v1].add(v2)
        vertex_neighbors[v2].add(v1)

    if not vertex_neighbors:
        return "other-open"

    vertex_degree = {v: len(neigh) for v, neigh in vertex_neighbors.items()}
    open_vertices = [v for v, deg in vertex_degree.items() if deg == 1]

    if not open_vertices:
        return "closed"

    all_on_axis = True
    for v in open_vertices:
        x, y, z = geompy.PointCoordinates(v)
        if abs(x) > axis_tol:
            all_on_axis = False
            break

    return "axis-open" if all_on_axis else "other-open"
def match_component_to_wire_by_vertices(comp_vertices,
                                        comp_type,
                                        wire_cache,
                                        dist_tol=1e-3,
                                        min_frac=0.3):
    """
    Match a BC component to a domain wire using vertex proximity,
    restricted to wires with the same open/closed/axis-open type.
    """
    if not comp_vertices:
        return None, None, 0, None

    # BC type filter: only consider wires of the same type
    allowed_types = []
    if comp_type == "closed":
        allowed_types = ["closed"]
    elif comp_type == "axis-open":
        allowed_types = ["axis-open"]
    else:
        # for now: allow all, but you can choose to drop 'other-open' entirely
        allowed_types = ["closed", "axis-open", "other-open"]

    bc_coords = [geompy.PointCoordinates(v) for v in comp_vertices]

    best_wire   = None
    best_index  = None
    best_score  = 0
    best_avg    = None

    for entry in wire_cache:
        if entry["type"] not in allowed_types:
            continue

        w     = entry["wire"]
        idx   = entry["index"]
        vlist = entry["vertices"]

        matches = 0
        dsum    = 0.0

        for (cx, cy, cz) in bc_coords:
            mind = None
            for v_wire, (wx, wy, wz) in vlist:
                dx = cx - wx
                dy = cy - wy
                dz = cz - wz
                d  = math.sqrt(dx*dx + dy*dy + dz*dz)
                if mind is None or d < mind:
                    mind = d
            if mind is not None and mind < dist_tol:
                matches += 1
                dsum    += mind

        if matches == 0:
            continue

        frac = matches / float(len(bc_coords))
        if frac < min_frac:
            continue

        avg = dsum / float(matches)

        if (matches > best_score) or (matches == best_score and
                                      (best_avg is None or avg < best_avg)):
            best_score = matches
            best_avg   = avg
            best_wire  = w
            best_index = idx

    return best_wire, best_index, best_score, best_avg
def make_wire_from_edges(edge_list, name=None):
    """Build a wire from a list of EDGE objects."""
    try:
        if not edge_list:
            raise RuntimeError(f"No edges provided to build wire {name or ''}")
        w = geompy.MakeWire(edge_list)
        if name:
            geompy.addToStudy(w, name)
        return w
    except:
        return geompy.MakeCompound(edge_list)

# ------------------------------------------------------------
# Normalization: always end up with a list of EDGE objects
# ------------------------------------------------------------
def normalize_to_edges(obj, context_name):
    """
    Ensure we always have a list of EDGE objects.

    obj can be:
      - list/tuple of EDGE objects
      - a GEOM shape that contains EDGE sub-shapes (WIRE / COMPOUND / FACE / etc)
    """
    if isinstance(obj, (list, tuple)):
        edges = list(obj)
    else:
        # Single GEOM object: extract edges from it
        edges = geompy.ExtractShapes(obj, geompy.ShapeType["EDGE"], True)

    if not edges:
        print(f"[normalize_to_edges] WARNING: no edges found in {context_name}")
    else:
        print(f"[normalize_to_edges] {context_name}: {len(edges)} EDGE(s)")

    return edges


# -----------------------------------------------
# Low-level: build a wire from a set of edges
# -----------------------------------------------
def build_wire_from_edges(edges_or_shape, wire_name, auto_connect=True):
    """
    Build a GEOM wire from an input set of edges.

    edges_or_shape can be:
      - list/tuple of EDGE objects
      - a GEOM shape that contains EDGE subshapes (COMPOUND/WIRE/FACE/etc)

    Returns:
        wire (GEOM_Object) whose edges you can inspect with ExtractShapes.
    """
    # Normalise to a list of EDGE objects
    if isinstance(edges_or_shape, (list, tuple)):
        edges = list(edges_or_shape)
    else:
        edges = geompy.ExtractShapes(edges_or_shape, geompy.ShapeType["EDGE"], True)

    if not edges:
        raise RuntimeError(f"build_wire_from_edges: no edges for {wire_name}")

    # First attempt: direct MakeWire(list_of_edges)
    try:
        w = geompy.MakeWire(edges)
        geompy.addToStudy(w, wire_name)
        return w
    except Exception as e1:
        if not auto_connect:
            raise

        # Fallback: let SALOME reconnect a compound of all edges
        comp = geompy.MakeCompound(edges)
        try:
            w = geompy.MakeWire(comp, True)  # True = auto reconnect
            geompy.addToStudy(w, wire_name + "_auto")
            return w
        except Exception as e2:
            print(f"build_wire_from_edges: failed for {wire_name}: {e1} / {e2}")
            raise


# -----------------------------------------------
# Mid-level: make a face with holes from wires
# -----------------------------------------------
def make_face_with_holes(
    outer_edges,
    list_of_hole_edge_lists,
    face_name,
    publish_face_edges=True,
):
    """
    Build a planar face from:
      - outer_edges: edges or a shape for the outer loop
      - list_of_hole_edge_lists: iterable of edge sets (or shapes), each one loop

    Only the resulting face and its own edges matter here.
    """

    # Outer wire
    outer_wire = build_wire_from_edges(
        edges_or_shape=outer_edges,
        wire_name=f"{face_name}_outer_wire",
        auto_connect=True,
    )

    # Hole wires
    hole_wires = []
    for i, hole_obj in enumerate(list_of_hole_edge_lists, start=1):
        if not hole_obj:
            continue
        try:
            hw = build_wire_from_edges(
                edges_or_shape=hole_obj,
                wire_name=f"{face_name}_hole_{i}_wire",
                auto_connect=True,
            )
            hole_wires.append(hw)
        except Exception as e:
            print(
                f"make_face_with_holes: skipping hole loop {i} for {face_name} "
                f"(wire build failed: {e})"
            )

    all_wires = [outer_wire] + hole_wires

    # Final face
    face = geompy.MakeFaceWires(all_wires, 1)  # 1 = planar
    geompy.addToStudy(face, face_name)

    # Optionally expose the *face’s* own edges for inspection
    if publish_face_edges:
        face_edges = geompy.ExtractShapes(face, geompy.ShapeType["EDGE"], True)
        for i, e in enumerate(face_edges, start=1):
            geompy.addToStudy(e, f"{face_name}_edge_{i}")

    return face




if __name__ == '__main__':
    model.begin()

    partSet = model.moduleDocument()

    makeface = True
    debug = False

    # Make the cryostat
    Cryostat = model.addPart(partSet)
    Cryostat_doc = Cryostat.document()
    p = GeometryParams()

    # If you want to "remove" a component, set it to None instead of calling build_...
    LXeCryostat = build_lower_cryostat(Cryostat_doc, p, makeface)
    GXeCryostat = build_upper_cryostat(Cryostat_doc, p, makeface)


    # Any of the below can be set to none to debug what is breaking the geometry
    # If all electrodes are none it wont work tho, each geometry needs at least one object
    # intersecting it for cutting

    # Electrodes
    FieldCage      = None if debug else build_field_cage(Cryostat_doc, p, makeface)
    FieldCageGuard = None if debug else build_field_cage_guard(Cryostat_doc, p, makeface)
    Bell           = build_bell(Cryostat_doc, p, makeface)
    Gate           = None if debug else build_gate(Cryostat_doc, p, makeface)
    Anode          = None if debug else build_anode(Cryostat_doc, p, makeface)
    Cathode        = build_cathode(Cryostat_doc, p, makeface)  # TODO First wire in wrong spot i think
    TopScreen      = None if debug else build_topScreen(Cryostat_doc, p, makeface)
    BottomScreen   = None if debug else build_bottomScreen(Cryostat_doc, p, makeface)  # TODO Needs the pin
    CopperRing     = None if debug else build_copper_ring(Cryostat_doc, p, makeface)
    # Too many vertices for the kernel so had to seperate them out 
    BottomPMTS     = None if debug else build_bottom_pmts(Cryostat_doc, p, makeface)
    TopPMTS        = None if debug else build_top_pmts(Cryostat_doc, p, makeface)

    ## PTFE Components
    ## TODO Missing the ptfe between copper ring and gate
    GateInsulatingFrame      = build_gate_insulating_frame(Cryostat_doc, p, makeface)
    AnodeInsulatingFrame     = build_anode_insulating_frame(Cryostat_doc, p, makeface)
    CathodeInsulatingFrame   = build_cathode_insulating_frame(Cryostat_doc, p, makeface)
    TopScreenInsulatingFrame = build_topscreen_insulating_frame(Cryostat_doc, p, makeface)
    BottomStackInsulation    = build_bottomstack_insulating_frame(Cryostat_doc, p, makeface)
    BottomPMTReflectors      = build_bottom_PMT_reflectors(Cryostat_doc, p, makeface)
    TopPMTReflectors         = build_top_PMT_reflectors(Cryostat_doc, p, makeface)
    wallPTFE                 = BuildPTFEWall(Cryostat_doc, p, makeface)

    model.do()

    # Helper: safe face selections (works even if feature is None)
    def safe_face_selections(feature):
        return face_feature_to_selections(feature) if feature is not None else []

    # ------------------------ We make face selections for each group ----------------------
    selec_Bell           = safe_face_selections(Bell)
    selec_Gate           = safe_face_selections(Gate)
    selec_Anode          = safe_face_selections(Anode)
    selec_Cathode        = safe_face_selections(Cathode)
    selec_TopScreen      = safe_face_selections(TopScreen)
    selec_BottomScreen   = safe_face_selections(BottomScreen)
    selec_CopperRing     = safe_face_selections(CopperRing)
    selec_FieldCage      = safe_face_selections(FieldCage)
    selec_FieldCageGuard = safe_face_selections(FieldCageGuard)
    selec_BottomPMTS     = safe_face_selections(BottomPMTS)
    selec_TopPMTS        = safe_face_selections(TopPMTS)

    # Per-PTFE selection lists
    selec_GateInsulatingFrame      = safe_face_selections(GateInsulatingFrame)
    selec_AnodeInsulatingFrame     = safe_face_selections(AnodeInsulatingFrame)
    selec_CathodeInsulatingFrame   = safe_face_selections(CathodeInsulatingFrame)
    selec_TopScreenInsulatingFrame = safe_face_selections(TopScreenInsulatingFrame)
    selec_BottomStackInsulation    = safe_face_selections(BottomStackInsulation)
    selec_BottomPMTReflectors      = safe_face_selections(BottomPMTReflectors)
    selec_TopPMTReflectors         = safe_face_selections(TopPMTReflectors)
    selec_wallPTFE                 = safe_face_selections(wallPTFE)

    # Background volumes (these are normally required; if you set one to None, this will just be [])
    selec_GXeCryostat = safe_face_selections(GXeCryostat)
    selec_LXeCryostat = safe_face_selections(LXeCryostat)

     # ----------------------------------------------------------------------
    # Per-object XAO exports (raw faces; no Shaper cut/boolean)
    # ----------------------------------------------------------------------

    def make_face_group(doc, selections, base_name):
        """Create a FACE group from a list of selections and wrap it in a GroupShape."""
        if not selections:
            return None, None

        g = model.addGroup(doc, "FACE", selections)
        g.setName(f"{base_name}Region")
        g.result().setName(f"{base_name}Region")

        # Wrap in a COMPOUND GroupShape (consistent with your PTFE pattern)
        gs = model.addGroupShape(doc, [model.selection("COMPOUND", f"{base_name}Region")])
        gs.setName(base_name)
        gs.result().setName(base_name)

        return g, gs
    
    # Cryostat envelopes
    GXe_region, GXe_group = make_face_group(Cryostat_doc, selec_GXeCryostat, "GXeCryostat")
    LXe_region, LXe_group = make_face_group(Cryostat_doc, selec_LXeCryostat, "LXeCryostat")

    # Electrodes (single logical object each)
    Bell_region,         Bell_group         = make_face_group(Cryostat_doc, selec_Bell,         "Bell")
    Gate_region,         Gate_group         = make_face_group(Cryostat_doc, selec_Gate,         "Gate")
    Anode_region,        Anode_group        = make_face_group(Cryostat_doc, selec_Anode,        "Anode")
    Cathode_region,      Cathode_group      = make_face_group(Cryostat_doc, selec_Cathode,      "Cathode")
    TopScreen_region,    TopScreen_group    = make_face_group(Cryostat_doc, selec_TopScreen,    "TopScreen")
    BottomScreen_region, BottomScreen_group = make_face_group(Cryostat_doc, selec_BottomScreen, "BottomScreen")
    CopperRing_region,   CopperRing_group   = make_face_group(Cryostat_doc, selec_CopperRing,   "CopperRing")
    BottomPMTs_region,   BottomPMTs_group   = make_face_group(Cryostat_doc, selec_BottomPMTS,   "BottomPMTs")
    TopPMTs_region,      TopPMTs_group      = make_face_group(Cryostat_doc, selec_TopPMTS,      "TopPMTs")

    # PTFE parts
    GateIns_region,      GateIns_group      = make_face_group(Cryostat_doc, selec_GateInsulatingFrame,      "GateInsulatingFrame")
    AnodeIns_region,     AnodeIns_group     = make_face_group(Cryostat_doc, selec_AnodeInsulatingFrame,     "AnodeInsulatingFrame")
    CathodeIns_region,   CathodeIns_group   = make_face_group(Cryostat_doc, selec_CathodeInsulatingFrame,   "CathodeInsulatingFrame")
    TopScreenIns_region, TopScreenIns_group = make_face_group(Cryostat_doc, selec_TopScreenInsulatingFrame, "TopScreenInsulatingFrame")
    BottomStackIns_region, BottomStackIns_group = make_face_group(Cryostat_doc, selec_BottomStackInsulation, "BottomStackInsulatingFrame")
    BottomPMTRefl_region, BottomPMTRefl_group = make_face_group(Cryostat_doc, selec_BottomPMTReflectors,    "BottomPMTReflectors")
    TopPMTRefl_region,    TopPMTRefl_group    = make_face_group(Cryostat_doc, selec_TopPMTReflectors,       "TopPMTReflectors")
    wallPTFE_region,      wallPTFE_group      = make_face_group(Cryostat_doc, selec_wallPTFE,               "wallPTFE")

    # ----------------------------------------------------------------------
    # FieldCage / FieldCageGuard: one group per face, with index suffix
    # ----------------------------------------------------------------------
    FieldCage_groups = []
    for i, sel in enumerate(selec_FieldCage, start=1):
        if not sel:
            continue
        base = f"FieldCage_{i}"
        g = model.addGroup(Cryostat_doc, "FACE", [sel])
        g.setName(f"{base}Region")
        g.result().setName(f"{base}Region")

        gs = model.addGroupShape(
            Cryostat_doc,
            [model.selection("COMPOUND", f"{base}Region")]
        )
        gs.setName(base)
        gs.result().setName(base)

        FieldCage_groups.append(gs)

    FieldCageGuard_groups = []
    for i, sel in enumerate(selec_FieldCageGuard, start=1):
        if not sel:
            continue
        base = f"FieldCageGuard_{i}"
        g = model.addGroup(Cryostat_doc, "FACE", [sel])
        g.setName(f"{base}Region")
        g.result().setName(f"{base}Region")

        gs = model.addGroupShape(
            Cryostat_doc,
            [model.selection("COMPOUND", f"{base}Region")]
        )
        gs.setName(base)
        gs.result().setName(base)

        FieldCageGuard_groups.append(gs)

    model.do()

    def xao_path_for(name):
        return f"/tmp/{name}.xao"

    xao_files = {}  # logical name -> path

    # Simple single objects: export by GroupShape name
    group_shapes = [
        GXe_group, LXe_group,
        Bell_group, Gate_group, Anode_group, Cathode_group,
        TopScreen_group, BottomScreen_group, CopperRing_group,
        BottomPMTs_group, TopPMTs_group,
        GateIns_group, AnodeIns_group, CathodeIns_group,
        TopScreenIns_group, BottomStackIns_group,
        BottomPMTRefl_group, TopPMTRefl_group, wallPTFE_group,
    ]

    for gs in group_shapes:
        if gs is None:
            continue
        base_name = gs.name()              # logical name: "Cathode", "LXeCryostat", ...
        res_name  = gs.result().name()     # actual COMPOUND result name to select
        path = xao_path_for(base_name)
        sel  = model.selection("COMPOUND", res_name)
        model.exportToXAO(Cryostat_doc, path, sel, "XAO")
        xao_files[base_name] = path
        if debug:
            print(f"Exported {base_name} (result '{res_name}') -> {path}")

    # FieldCage / FieldCageGuard per-ring exports
    for gs in FieldCage_groups + FieldCageGuard_groups:
        if gs is None:
            continue
        base_name = gs.name()              # e.g. "FieldCage_1"
        res_name  = gs.result().name()
        path = xao_path_for(base_name)
        sel  = model.selection("COMPOUND", res_name)
        model.exportToXAO(Cryostat_doc, path, sel, "XAO")
        xao_files[base_name] = path
        if debug:
            print(f"Exported {base_name} (result '{res_name}') -> {path}")

    model.end()
    print("Finished SHAPER Component (primitive/group exports)")
    # FIXME 

    # ------------ We export to the old geometry module because i cant figure it out ------------
    # First export to geom make one parent module and then make everything a subgroup so we 
    # have volume and edge elements 
    
    ###
    ### GEOM component
    ###

    import GEOM
    from salome.geom import geomBuilder
    import SALOMEDS

    geompy = geomBuilder.New()

    # Reference axes (optional)
    O  = geompy.MakeVertex(0, 0, 0)
    OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
    OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

    geompy.addToStudy(O,  'O')
    geompy.addToStudy(OX, 'OX')
    geompy.addToStudy(OY, 'OY')
    geompy.addToStudy(OZ, 'OZ')

    # ----------------------------------------------------------------------
    # Import all XAO primitives exported from Shaper
    # xao_files: dict from Shaper stage, name -> path (e.g. "Cathode": "/tmp/Cathode.xao")
    # For some reason not all objects are faces - we simply import edges and make all teh faces again 
    # ----------------------------------------------------------------------

    imported_edges  = {}  # logical name -> [EDGE, EDGE, ...]
    imported_raw    = {}  # logical name -> raw imported GEOM shape (for inspection only)

    for name, path in xao_files.items():
        imported, shape, sub, groups, fields = geompy.ImportXAO(path)
        imported_raw[name] = shape

        # Extract all edges from this XAO
        edges = geompy.ExtractShapes(shape, geompy.ShapeType["EDGE"], True)
        if not edges:
            print(f"Warning: no edges found in XAO for '{name}'")
            imported_edges[name] = []
        else:
            imported_edges[name] = edges

        # Optional: publish raw imported shape for debugging
        geompy.addToStudy(shape, f"RAW_{name}")

    print("GEOM: collected edges for:", list(imported_edges.keys()))

    
    # ----------------------------------------------------------------------
    # 3) Classify edge-sets into LXe / GXe / PTFE
    #    (YOU must fill these lists according to where each part lives.)
    # ----------------------------------------------------------------------

    # Names for the cryostat outlines from Shaper
    LXe_outline_name = "LXeCryostat"
    GXe_outline_name = "GXeCryostat"

    # PTFE parts (all holes inside both LXe/GXe, depending on geometry)
    ptfe_inLXe = [
        "GateInsulatingFrame",
        "CathodeInsulatingFrame",
        "BottomStackInsulatingFrame",
        "BottomPMTReflectors",
        "wallPTFE",
    ]
    ptfe_inGXe = [
        "CathodeInsulatingFrame",
        "TopPMTReflectors",
    ]

    # Electrodes that sit in LXe (holes in LXeRegion)
    electrodes_in_LXe = [
        "Gate",
        "Cathode",
        "BottomPMTs",
        "BottomScreen",
        "CopperRing",
        # ... add more if needed
    ]

    # Electrodes that sit in GXe (holes in GXeRegion)
    electrodes_in_GXe = [
        "Anode",
        "TopScreen",
        "TopPMTs",
        "Bell",
        # ... etc.
    ]
    # ---------------------- LXeRegion ----------------------------
    outer_LXe_edges = imported_edges.get(LXe_outline_name, [])

    LXe_hole_edge_lists = []
    for n in ptfe_inLXe + electrodes_in_LXe:
        if n not in imported_edges:
            continue
        LXe_hole_edge_lists.append(imported_edges[n])

    LXeRegion = make_face_with_holes(
        outer_edges      = outer_LXe_edges,          # what you used before
        list_of_hole_edge_lists = LXe_hole_edge_lists,    # list of loops per PTFE/electrode
        face_name        = "LXeRegion",
        publish_face_edges = True,                  # handy while debugging
    )




    # ---------------------- GXeRegion ----------------------------
    outer_GXe_edges = imported_edges.get(GXe_outline_name, [])

    GXe_hole_edge_lists = []
    for n in ptfe_inGXe + electrodes_in_GXe:
        if n not in imported_edges:
            continue
        GXe_hole_edge_lists.append(imported_edges[n])

    GXeRegion = make_face_with_holes(
        outer_edges      = outer_GXe_edges,
        list_of_hole_edge_lists = GXe_hole_edge_lists,
        face_name        = "GXeRegion",
        publish_face_edges = True,
    )
    raise Exception()

    # ---------------------- PTFERegion ---------------------------
    # You can either:
    #  - build separate PTFE faces per object, then MakeCompound
    #  - or a single PTFE face if topology is simple.

    ptfe_faces = []
    for n in ptfe_inGXe + ptfe_inLXe:
        edges = imported_edges.get(n, [])
        if not edges:
            continue
        w = make_wire_from_edges(edges, name=f"{n}_wire")
        f = geompy.MakeFaceWires([w], 1)
        geompy.addToStudy(f, f"{n}_face")
        ptfe_faces.append(f)

    PTFERegion = geompy.MakeCompound(ptfe_faces) if ptfe_faces else None
    if PTFERegion:
        geompy.addToStudy(PTFERegion, "PTFERegion")


    # ----------------------------------------------------------------------
    # Build BC edge geometries directly from imported_edges
    # ----------------------------------------------------------------------

    bc_geom = {}  # BC_<ElectrodeName> -> GEOM edge compound
    electrode_names = electrodes_in_GXe + electrodes_in_LXe

    for elec_name in electrode_names:
        edges = imported_edges.get(elec_name, [])
        if not edges:
            print(f"Warning: no edges stored for electrode '{elec_name}'")
            continue

        # One compound of all electrode edges for this BC
        bc_shape = geompy.MakeCompound(edges)

        bc_name = f"BC_{elec_name}"       # e.g. BC_Cathode, BC_FieldCage_3
        bc_geom[bc_name] = bc_shape

        # Optional: publish for inspection in the GEOM viewer
        geompy.addToStudy(bc_shape, bc_name + "_Geom")

    # Meshing domain: materials only
    meshing_domain_shapes = [s for s in (LXeRegion, GXeRegion, PTFERegion) if s is not None]
    MESHING_group = geompy.MakeCompound(meshing_domain_shapes)
    geompy.addToStudy(MESHING_group, "MESHING_group")
    

     # Parent object for meshing: ONLY fluids + PTFE
    meshing_domain_shapes = [s for s in (LXeRegion, GXeRegion, PTFERegion) if s is not None]
    MESHING_group = geompy.MakeCompound(meshing_domain_shapes)
    geompy.addToStudy(MESHING_group, "MESHING_group")

    print("GEOM: MESHING_group (LXe + GXe + PTFE) constructed.")

    # ------------------------------------------------------------------
    # Clip BC raw edge geometries to the meshing domain using MakeSection
    #
    # bc_geom maps:
    #   bc_name -> raw electrode edge compound (may extend outside domain)
    #
    # We replace each entry by the intersection curves between that
    # BC geometry and the meshing domain.
    # ------------------------------------------------------------------

    bc_geom_clipped = {}

    for bc_name, raw_bc_shape in bc_geom.items():
        try:
            # Intersection of BC geometry with the meshing domain
            sec = geompy.MakeSection(MESHING_group, raw_bc_shape)
        except Exception as e:
            print(f"Warning: MakeSection failed for BC '{bc_name}': {e}")
            continue

        # Extract intersection edges
        sec_edges = geompy.ExtractShapes(sec, geompy.ShapeType["EDGE"], True)

        if not sec_edges:
            print(f"Warning: BC '{bc_name}' has no intersection with the mesh domain")
            continue

        # Repack as a compound of boundary edges
        clipped = geompy.MakeCompound(sec_edges)
        bc_geom_clipped[bc_name] = clipped

        # Optional: publish for visual inspection
        geompy.addToStudy(clipped, bc_name + "_Clipped")

    # Replace the original dict
    bc_geom = bc_geom_clipped

    print("GEOM: BC edge shapes clipped to MESHING_group.")


    ###
    ### SMESH component
    ###

    import SMESH, SALOMEDS
    from salome.smesh import smeshBuilder

    smesh = smeshBuilder.New()

    Mesh_1 = smesh.Mesh(MESHING_group)

    NETGEN_2D = Mesh_1.Triangle(smeshBuilder.NETGEN_1D2D)
    hyp2d = NETGEN_2D.Parameters()
    hyp2d.SetSecondOrder(0)
    hyp2d.SetFineness(1)            # medium
    hyp2d.SetOptimize(1)
    hyp2d.SetUseSurfaceCurvature(1)

    # ------------------ Region face groups -------------------
    if GXeRegion is not None:
        GXeRegion_1 = Mesh_1.GroupOnGeom(GXeRegion, 'GXeRegion', SMESH.FACE)
        smesh.SetName(GXeRegion_1, 'GXeRegion')

    if PTFERegion is not None:
        PTFERegion_1 = Mesh_1.GroupOnGeom(PTFERegion, 'PTFERegion', SMESH.FACE)
        smesh.SetName(PTFERegion_1, 'PTFERegion')

    if LXeRegion is not None:
        LXeRegion_1 = Mesh_1.GroupOnGeom(LXeRegion, 'LXeRegion', SMESH.FACE)
        smesh.SetName(LXeRegion_1, 'LXeRegion')

    smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 1D-2D')
    smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

    # ------------------ BC edge groups -----------------------
    bc_mesh_groups = {}

    for bc_name, bc_shape in bc_geom.items():
        # All mesh edges lying on those geometry edges
        grp = Mesh_1.GroupOnGeom(bc_shape, bc_name, SMESH.EDGE)
        smesh.SetName(grp, bc_name)
        bc_mesh_groups[bc_name] = grp

    # Now compute and export
    #Mesh_1.Compute()
    Mesh_1.ExportMED(
        "/home/felix/MFEMElectrostatics/geometries/SR3nTCOMSOLVerification/mesh.med",
        auto_groups=False,
        minor=42,
    )
    print("Finished Meshing")
    if salome.sg.hasDesktop():
        salome.sg.updateObjBrowser()

    """    Convert to gmsh 2.2 and dump conversion 
    gmsh -0 mesh.med -format msh2 -o mesh22.msh
    awk '/\$PhysicalNames/{flag=1;next}/\$EndPhysicalNames/{flag=0}flag' mesh22.msh \
      > phys_map.txt
    python3 make_elements.py
    # And then merge config and config_autogen.py
    """
