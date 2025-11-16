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
    FieldCage      = None#build_field_cage(Cryostat_doc, p, makeface)
    FieldCageGuard = None#build_field_cage_guard(Cryostat_doc, p, makeface)
    Bell           = build_bell(Cryostat_doc, p, makeface)
    Gate           = None#build_gate(Cryostat_doc, p, makeface)
    Anode          = None#build_anode(Cryostat_doc, p, makeface)
    Cathode        = build_cathode(Cryostat_doc, p, makeface)  # TODO First wire in wrong spot i think
    TopScreen      = None#build_topScreen(Cryostat_doc, p, makeface)
    BottomScreen   = None#build_bottomScreen(Cryostat_doc, p, makeface)  # TODO Needs the pin
    CopperRing     = None#build_copper_ring(Cryostat_doc, p, makeface)
    # Too many vertices for the kernel so had to seperate them out 
    BottomPMTS     = None#build_bottom_pmts(Cryostat_doc, p, makeface)
    TopPMTS        = None#build_top_pmts(Cryostat_doc, p, makeface)

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
    # Generate tools lists 
    # ----------------------------------------------------------------------
    all_electrode = (
        selec_Bell
        + selec_Gate
        + selec_Anode
        + selec_Cathode
        + selec_TopScreen
        + selec_BottomScreen
        + selec_CopperRing
        + selec_FieldCage
        + selec_FieldCageGuard
        + selec_BottomPMTS
        + selec_TopPMTS
    )

    all_ptfe = (
        selec_GateInsulatingFrame
        + selec_AnodeInsulatingFrame
        + selec_CathodeInsulatingFrame
        + selec_TopScreenInsulatingFrame
        + selec_BottomStackInsulation
        + selec_BottomPMTReflectors
        + selec_TopPMTReflectors
        + selec_wallPTFE
    )

    # ----------------------------------------------------------------------
    # Cut objects out of background using only selections
    # Also make a union of PTFE 
    # ----------------------------------------------------------------------
    main_sel = selec_GXeCryostat + selec_LXeCryostat
    tool_sels = all_electrode + all_ptfe

    cutResult = model.addCut(
        Cryostat_doc,
        main_sel,   # list[ModelHighAPI_Selection]
        tool_sels,  # list[ModelHighAPI_Selection]
        keepSubResults=True
    )

    # -------------------------------     Set Material Group Names --------------------------------
    GXe_group = model.addGroup(
        Cryostat_doc,
        "FACE",
        [model.selection("FACE", "Cut_1_1")]
    )
    GXe_group.setName("GXeGroup")
    GXe_group.result().setName("GXeGroup")

    LXe_group = model.addGroup(
        Cryostat_doc,
        "FACE",
        [model.selection("FACE", "Cut_1_2")]
    )
    LXe_group.setName("LXeGroup")
    LXe_group.result().setName("LXeGroup")

    GXeVol, LXeVol = cutResult.results()
    cutResult.setName("VolumeCut")
    LXeVol.setName("LXeRegion")
    GXeVol.setName("GXeRegion")

    # Create PTFE region group only if there is something in all_ptfe
    if all_ptfe:
        PTFE_region = model.addGroup(
            Cryostat_doc,
            "FACE",
            all_ptfe
        )
        PTFE_region.setName("PTFERegion")
        PTFE_region.result().setName("PTFERegion")

        PTFE_group = model.addGroupShape(
            Cryostat_doc, [model.selection("COMPOUND", "PTFERegion")]
        )
        PTFE_group.setName("PTFE_Group")
        PTFE_group.result().setName("PTFE_Group")
    else:
        PTFE_region = None
        PTFE_group = None

    model.do()

    # ---------------------------------- Set BC group Names ------------------------------------------

    def add_edge_group(doc, selections, name):
        """Create an EDGE group only if there are selections; otherwise return None."""
        if not selections:
            return None
        g = model.addGroup(doc, "EDGE", selections)
        g.setName(name)
        g.result().setName(name)
        return g

    BellBC         = add_edge_group(Cryostat_doc, selec_Bell,         "BC_Bell")
    GateBC         = add_edge_group(Cryostat_doc, selec_Gate,         "BC_Gate")
    AnodeBC        = add_edge_group(Cryostat_doc, selec_Anode,        "BC_Anode")
    CathodeBC      = add_edge_group(Cryostat_doc, selec_Cathode,      "BC_Cathode")
    TopScreenBC    = add_edge_group(Cryostat_doc, selec_TopScreen,    "BC_TopScreen")
    BottomScreenBC = add_edge_group(Cryostat_doc, selec_BottomScreen, "BC_BottomScreen")
    CopperRingBC   = add_edge_group(Cryostat_doc, selec_CopperRing,   "BC_CopperRing")
    PMTBC          = add_edge_group(Cryostat_doc, selec_BottomPMTS + selec_TopPMTS, "BC_PMTs")

    FieldCageBC = []
    for i, sel in enumerate(selec_FieldCage, start=1):
        name = f"BC_FieldCage_Ring_{i}"
        grp = add_edge_group(Cryostat_doc, [sel], name)
        if grp is not None:
            FieldCageBC.append(grp)

    FieldCageGuardBC = []
    for i, sel in enumerate(selec_FieldCageGuard, start=1):
        name = f"BC_FieldCageGuard_Ring_{i}"
        grp = add_edge_group(Cryostat_doc, [sel], name)
        if grp is not None:
            FieldCageGuardBC.append(grp)

    model.do()

    ### Create Exports since apparently it has to go by file
    # Volumes
    Export_1 = model.exportToXAO(
        Cryostat_doc,
        '/tmp/shaper_6c6wg3w7.xao',
        model.selection("COMPOUND", "GXeRegion"),
        'XAO'
    )
    Export_2 = model.exportToXAO(
        Cryostat_doc,
        '/tmp/shaper_b_d4_v11.xao',
        model.selection("COMPOUND", "LXeRegion"),
        'XAO'
    )

    if PTFE_group is not None:
        Export_3 = model.exportToXAO(
            Cryostat_doc,
            '/tmp/shaper_fyn9_2he.xao',
            model.selection("COMPOUND", "PTFE_Group"),
            'XAO'
        )
    else:
        Export_3 = None

    # Boundaries 
    def tmp_xao_path(prefix="shaper"):
        return f"/tmp/{prefix}_{uuid.uuid4().hex[:8]}.xao"

    all_bc_groups = [
        BellBC, GateBC, AnodeBC, CathodeBC, TopScreenBC, BottomScreenBC,
        CopperRingBC, PMTBC, *FieldCageBC, *FieldCageGuardBC
    ]
    # Drop None groups (for disabled/empty components)
    all_bc_groups = [g for g in all_bc_groups if g is not None]

    bc_xao_files = []
    for grp in all_bc_groups:
        xao_path = tmp_xao_path(grp.name())
        sel = model.selection("COMPOUND", grp.name())
        model.exportToXAO(Cryostat_doc, xao_path, sel, "XAO")
        bc_xao_files.append((grp.name(), xao_path))

    model.end()

    # ------------ We export to the old geometry module because i cant figure it out ------------
    # First export to geom make one parent module and then make everything a subgroup so we 
    # have volume and edge elements 
    
    ###
    ### GEOM component
    ###

    import GEOM
    from salome.geom import geomBuilder
    import math
    import SALOMEDS


    geompy = geomBuilder.New()

    O = geompy.MakeVertex(0, 0, 0)
    OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
    OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
    # Volumes
    (imported, GXeRegion, [], [GXeGroup], []) = geompy.ImportXAO("/tmp/shaper_6c6wg3w7.xao")
    (imported, LXeRegion, [], [LXeGroup], []) = geompy.ImportXAO("/tmp/shaper_b_d4_v11.xao")
    (imported, PTFE_Group, [], [], []) = geompy.ImportXAO("/tmp/shaper_fyn9_2he.xao")
    faces_GXe = geompy.MakeCompound(geompy.ExtractShapes(GXeRegion, geompy.ShapeType["FACE"], True))
    faces_LXe = geompy.MakeCompound(geompy.ExtractShapes(LXeRegion, geompy.ShapeType["FACE"], True))
    faces_PTFE = geompy.MakeCompound(geompy.ExtractShapes(PTFE_Group, geompy.ShapeType["FACE"], True))
    # BC's 
    bc_geom = {}
    bc_edges = []  # Stores all extracted EDGE objects, purely for optional diagnostics
    for bc_name, xao_file in bc_xao_files:
        imported, bc_shape, sub_shapes, groups, fields = geompy.ImportXAO(xao_file)
        # Extract edges associated with this BC
        edges = geompy.ExtractShapes(bc_shape, geompy.ShapeType["EDGE"], True)
        if not edges:
            print(f"Warning: no edges extracted for BC '{bc_name}' from {xao_file}")
            continue
        # If multiple edges exist, make a compound; otherwise, use the single edge
        bc_edge_obj = edges[0] if len(edges) == 1 else geompy.MakeCompound(edges)
        # Add to the study with a meaningful name
        geompy.addToStudy(bc_edge_obj, bc_name + "_Geom")
        # Store for later: (geometry, optional XAO group)
        bc_geom[bc_name] = (bc_edge_obj, groups)
        # Accumulate edges for diagnostics/debug
        bc_edges.extend(edges)

    # Make parent object
    MESHING_group = geompy.MakeShell([GXeRegion, PTFE_Group, LXeRegion])
    #MESHING_group = geompy.MakeCompound([faces_GXe, faces_PTFE, faces_LXe])
    geompy.addToStudy( O, 'O' )
    geompy.addToStudy( OX, 'OX' )
    geompy.addToStudy( OY, 'OY' )
    geompy.addToStudy( OZ, 'OZ' )
    # Volumes
    geompy.addToStudy( GXeRegion, 'GXeRegion' )
    geompy.addToStudyInFather( GXeRegion, GXeGroup, 'GXeGroup' )
    geompy.addToStudy( LXeRegion, 'LXeRegion' )
    geompy.addToStudyInFather( LXeRegion, LXeGroup, 'LXeGroup' )
    geompy.addToStudy( PTFE_Group, 'PTFE_Group' )
    # Parent Object 
    geompy.addToStudy( MESHING_group, 'MESHING_group' )

    ###
    ### SMESH component
    ###

    import  SMESH, SALOMEDS
    from salome.smesh import smeshBuilder

    smesh = smeshBuilder.New()
    #smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                    # multiples meshes built in parallel, complex and numerous mesh edition (performance)
    Mesh_1 = smesh.Mesh(MESHING_group)
    algo1d = Mesh_1.Segment() 
    local_wire_length = min(p.TopScreenWireDiameter, p.GateWireDiameter, p.TopStackWireDiameter, p.CathodeWireDiameter, p.BottomScreenWireDiameter)
    local_wire_length *= math.pi *2 # Want 12 edges about the wire
    hyp1d = algo1d.LocalLength(local_wire_length)
    NETGEN_2D = Mesh_1.Triangle(smeshBuilder.NETGEN_2D)
    hyp2d  = NETGEN_2D.Parameters()
    #hyp2d.SetMaxSize(5.0 * local_wire_length)           # global max element size
    #hyp2d.SetMinSize(0.25 * local_wire_length)         # global min element size (0 disables)
    hyp2d.SetSecondOrder(0)
    hyp2d.SetFineness(1)  # medium granularity
    hyp2d.SetOptimize(1)
    hyp2d.SetUseSurfaceCurvature(1)

    #hyp2d.SetSecondOrder(0)        # 0 = first order, 1 = second order
    #hyp2d.SetOptimize(0)           # 1 = optimize, 0 = off
    #hyp2d.SetFineness(0)           # 0..4 (VeryCoarse..VeryFine), 2–3 is typical
    #hyp2d.SetUseSurfaceCurvature(1)  # refine near curvature
    #hyp2d.SetGrowthRate(0.3)       # 0..1, smaller -> smoother gradation
    #hyp2d.SetNbSegPerEdge(1)       # segments per edge (used with LocalLength)
    #hyp2d.SetNbSegPerRadius(2)     # segments per radius for curved entities
    #hyp2d.SetQuadAllowed(0)        # 0 = triangles only

    # ---------------------------------------------------------------------
    # Region groups: Faces
    # ---------------------------------------------------------------------
    GXeRegion_1 = Mesh_1.GroupOnGeom(GXeRegion, 'GXeRegion', SMESH.FACE)
    PTFE_Group_1 = Mesh_1.GroupOnGeom(PTFE_Group, 'PTFE_Group', SMESH.FACE)
    LXeRegion_1 = Mesh_1.GroupOnGeom(LXeRegion, 'LXeRegion', SMESH.FACE)

    smesh.SetName(LXeRegion_1, 'LXeRegion')
    smesh.SetName(PTFE_Group_1, 'PTFE_Group')
    smesh.SetName(GXeRegion_1, 'GXeRegion')

    #GXeRegion_2 = Mesh_1.GroupOnGeom(GXeRegion, 'GXeRegion', SMESH.NODE)
    #PTFE_Group_2 = Mesh_1.GroupOnGeom(PTFE_Group, 'PTFE_Group', SMESH.NODE)
    #LXeRegion_2 = Mesh_1.GroupOnGeom(LXeRegion, 'LXeRegion', SMESH.NODE)
    
    #smesh.SetName(GXeRegion_2, 'GXeRegion')
    #smesh.SetName(LXeRegion_2, 'LXeRegion')
    #smesh.SetName(PTFE_Group_2, 'PTFE_Group')

    smesh.SetName(algo1d.GetAlgorithm(),  'Algo 1D')
    smesh.SetName(hyp1d, "1D local length")
    smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 1D-2D')
    smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

    # ---------------------------------------------------------------------
    # Boundary-condition groups: edges on the meshed surface
    # ---------------------------------------------------------------------
    bc_mesh_groups = {}

    for bc_name, (bc_shape, groups) in bc_geom.items():
        # Change SMESH.EDGE to SMESH.FACE if your BC XAOs are faces
        grp = Mesh_1.GroupOnGeom(bc_shape, bc_name, SMESH.EDGE)
        #grp2 = Mesh_1.GroupOnGeom(bc_shape, bc_name, SMESH.NODE)
        smesh.SetName(grp, bc_name)
        #smesh.SetName(grp2, bc_name)
        bc_mesh_groups[bc_name] = grp


    #isDone = Mesh_1.Compute()
    #Mesh_1.ExportMED("/home/felix/MFEMElectrostatics/geometries/SR3nTCOMSOLVerification/mesh.med", 
    #                auto_groups=False, 
    #                minor=42, # MED version 4.2
    #                )

    if salome.sg.hasDesktop():
      salome.sg.updateObjBrowser()

    """    Convert to gmsh 2.2 and dump conversion 
    gmsh -0 mesh.med -format msh2 -o mesh22.msh
    awk '/\$PhysicalNames/{flag=1;next}/\$EndPhysicalNames/{flag=0}flag' mesh22.msh \
      > phys_map.txt
    python3 make_elements.py
    # And then merge config and config_autogen.py
    """
