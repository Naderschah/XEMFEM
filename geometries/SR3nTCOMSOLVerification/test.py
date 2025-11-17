#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/felix/MFEMElectrostatics/geometries/SR3nTCOMSOLVerification')

###
### SHAPER component
###

from SketchAPI import *

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

### Create Sketch
Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_1 = Sketch_1.addLine(0, -1.60824, 0, 0.004)
SketchLine_1.setName("AxisRect")
SketchLine_1.result().setName("SketchLine_1")

### Create SketchLine
SketchLine_2 = Sketch_1.addLine(0, 0.004, 0.73, 0.004)
SketchLine_2.setName("TopLine")
SketchLine_2.result().setName("SketchLine_2")

### Create SketchLine
SketchLine_3 = Sketch_1.addLine(0.73, 0.004, 0.73, -1.60824)
SketchLine_3.setName("RightWall")
SketchLine_3.result().setName("SketchLine_3")

### Create SketchArc
SketchArc_1 = Sketch_1.addArc(0, -0.7714007083794403, 7.347880794884119e-17, -1.97140070837944, 0.6244897959183675, -1.79610188179237, False)
SketchArc_1.setName("LowerCryostatArc")
SketchArc_1.result().setName("SketchArc_1")
SketchArc_1.results()[1].setName("SketchArc_1_2")

### Create SketchLine
SketchLine_4 = Sketch_1.addLine(0, -1.60824, 7.347880794884119e-17, -1.97140070837944)
SketchLine_4.setName("AxisClosingLine")
SketchLine_4.result().setName("SketchLine_4")

### Create SketchArc
SketchArc_2 = Sketch_1.addArc(0.51, -1.60824, 0.6244897959183675, -1.79610188179237, 0.73, -1.60824, False)
SketchArc_2.setName("UpperRoundingArc")
SketchArc_2.result().setName("SketchArc_2")
SketchArc_2.results()[1].setName("SketchArc_2_2")
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchLine_2.startPoint())
Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchLine_3.startPoint())
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchArc_2.endPoint())
Sketch_1.setCoincident(SketchArc_2.startPoint(), SketchArc_1.endPoint())
Sketch_1.setCoincident(SketchArc_1.startPoint(), SketchLine_4.endPoint())
Sketch_1.setCoincident(SketchLine_4.startPoint(), SketchLine_1.startPoint())
model.do()
Sketch_1.setName("LowerCryostatSketch")
Sketch_1.result().setName("LowerCryostatSketch")

### Create Face
Face_1 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "LowerCryostatSketch")])
Face_1.setName("LXeCryostat")
Face_1.result().setName("Face_1_1")

### Create Sketch
Sketch_2 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_5 = Sketch_2.addLine(0, 0.004, 0, 0.2588600000000001)

### Create SketchLine
SketchLine_6 = Sketch_2.addLine(0, 0.004, 0.73, 0.004)

### Create SketchLine
SketchLine_7 = Sketch_2.addLine(0.73, 0.004, 0.73, 0.2588600000000001)

### Create SketchArc
SketchArc_3 = Sketch_2.addArc(0, -0.5779792916205595, 0.6244897959183675, 0.4467218817923704, 7.347880794884119e-17, 0.6220207083794405, False)

### Create SketchArc
SketchArc_4 = Sketch_2.addArc(0.51, 0.2588600000000001, 0.73, 0.2588600000000001, 0.6244897959183675, 0.4467218817923704, False)

### Create SketchLine
SketchLine_8 = Sketch_2.addLine(0, 0.2588600000000001, 7.347880794884119e-17, 0.6220207083794405)
Sketch_2.setCoincident(SketchLine_5.endPoint(), SketchLine_8.startPoint())
Sketch_2.setCoincident(SketchLine_8.endPoint(), SketchArc_3.endPoint())
Sketch_2.setCoincident(SketchArc_3.startPoint(), SketchArc_4.endPoint())
Sketch_2.setCoincident(SketchArc_4.startPoint(), SketchLine_7.endPoint())
Sketch_2.setCoincident(SketchLine_7.startPoint(), SketchLine_6.endPoint())
Sketch_2.setCoincident(SketchLine_6.startPoint(), SketchLine_5.startPoint())
model.do()
Sketch_2.setName("UpperCryostatSketch")
Sketch_2.result().setName("UpperCryostatSketch")

### Create Face
Face_2 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "UpperCryostatSketch")])
Face_2.setName("GXeCryostat")
Face_2.result().setName("Face_2_1")

### Create Sketch
Sketch_3 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchCircle
SketchCircle_1 = Sketch_3.addCircle(0.668, -0.023, 0.001)

### Create SketchPoint
SketchPoint_1 = Sketch_3.addPoint(0.668, -0.023)

### Create SketchPoint
SketchPoint_2 = Sketch_3.addPoint(0.668, -0.034)

### Create SketchMultiTranslation
SketchMultiTranslation_1 = Sketch_3.addTranslation([SketchCircle_1.results()[1]], SketchPoint_1.coordinates(), SketchPoint_2.coordinates(), 4)
[SketchCircle_2, SketchCircle_3, SketchCircle_4] = SketchMultiTranslation_1.translatedList()

### Create SketchCircle
SketchCircle_5 = Sketch_3.addCircle(0.668, -0.067, 0.001)

### Create SketchPoint
SketchPoint_3 = Sketch_3.addPoint(0.668, -0.067)

### Create SketchPoint
SketchPoint_4 = Sketch_3.addPoint(0.668, -0.089)

### Create SketchMultiTranslation
SketchMultiTranslation_2 = Sketch_3.addTranslation([SketchCircle_5.results()[1]], SketchPoint_3.coordinates(), SketchPoint_4.coordinates(), 65)
[SketchCircle_6, SketchCircle_7, SketchCircle_8, SketchCircle_9, SketchCircle_10, SketchCircle_11, SketchCircle_12, SketchCircle_13, SketchCircle_14, SketchCircle_15, SketchCircle_16, SketchCircle_17, SketchCircle_18, SketchCircle_19, SketchCircle_20, SketchCircle_21, SketchCircle_22, SketchCircle_23, SketchCircle_24, SketchCircle_25, SketchCircle_26, SketchCircle_27, SketchCircle_28, SketchCircle_29, SketchCircle_30, SketchCircle_31, SketchCircle_32, SketchCircle_33, SketchCircle_34, SketchCircle_35, SketchCircle_36, SketchCircle_37, SketchCircle_38, SketchCircle_39, SketchCircle_40, SketchCircle_41, SketchCircle_42, SketchCircle_43, SketchCircle_44, SketchCircle_45, SketchCircle_46, SketchCircle_47, SketchCircle_48, SketchCircle_49, SketchCircle_50, SketchCircle_51, SketchCircle_52, SketchCircle_53, SketchCircle_54, SketchCircle_55, SketchCircle_56, SketchCircle_57, SketchCircle_58, SketchCircle_59, SketchCircle_60, SketchCircle_61, SketchCircle_62, SketchCircle_63, SketchCircle_64, SketchCircle_65, SketchCircle_66, SketchCircle_67, SketchCircle_68, SketchCircle_69] = SketchMultiTranslation_2.translatedList()

### Create SketchCircle
SketchCircle_70 = Sketch_3.addCircle(0.668, -1.486, 0.001)

### Create SketchPoint
SketchPoint_5 = Sketch_3.addPoint(0.668, -1.486)

### Create SketchPoint
SketchPoint_6 = Sketch_3.addPoint(0.668, -1.497)

### Create SketchMultiTranslation
SketchMultiTranslation_3 = Sketch_3.addTranslation([SketchCircle_70.results()[1]], SketchPoint_5.coordinates(), SketchPoint_6.coordinates(), 2)
[SketchCircle_71] = SketchMultiTranslation_3.translatedList()
model.do()
Sketch_3.setName("FieldCage")
Sketch_3.result().setName("FieldCage")

### Create Face
Face_3 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "FieldCage")])
Face_3.setName("FieldCage")
Face_3.result().setName("Face_3_1")
Face_3.results()[1].setName("Face_3_2")
Face_3.results()[2].setName("Face_3_3")
Face_3.results()[3].setName("Face_3_4")
Face_3.results()[4].setName("Face_3_5")
Face_3.results()[5].setName("Face_3_6")
Face_3.results()[6].setName("Face_3_7")
Face_3.results()[7].setName("Face_3_8")
Face_3.results()[8].setName("Face_3_9")
Face_3.results()[9].setName("Face_3_10")
Face_3.results()[10].setName("Face_3_11")
Face_3.results()[11].setName("Face_3_12")
Face_3.results()[12].setName("Face_3_13")
Face_3.results()[13].setName("Face_3_14")
Face_3.results()[14].setName("Face_3_15")
Face_3.results()[15].setName("Face_3_16")
Face_3.results()[16].setName("Face_3_17")
Face_3.results()[17].setName("Face_3_18")
Face_3.results()[18].setName("Face_3_19")
Face_3.results()[19].setName("Face_3_20")
Face_3.results()[20].setName("Face_3_21")
Face_3.results()[21].setName("Face_3_22")
Face_3.results()[22].setName("Face_3_23")
Face_3.results()[23].setName("Face_3_24")
Face_3.results()[24].setName("Face_3_25")
Face_3.results()[25].setName("Face_3_26")
Face_3.results()[26].setName("Face_3_27")
Face_3.results()[27].setName("Face_3_28")
Face_3.results()[28].setName("Face_3_29")
Face_3.results()[29].setName("Face_3_30")
Face_3.results()[30].setName("Face_3_31")
Face_3.results()[31].setName("Face_3_32")
Face_3.results()[32].setName("Face_3_33")
Face_3.results()[33].setName("Face_3_34")
Face_3.results()[34].setName("Face_3_35")
Face_3.results()[35].setName("Face_3_36")
Face_3.results()[36].setName("Face_3_37")
Face_3.results()[37].setName("Face_3_38")
Face_3.results()[38].setName("Face_3_39")
Face_3.results()[39].setName("Face_3_40")
Face_3.results()[40].setName("Face_3_41")
Face_3.results()[41].setName("Face_3_42")
Face_3.results()[42].setName("Face_3_43")
Face_3.results()[43].setName("Face_3_44")
Face_3.results()[44].setName("Face_3_45")
Face_3.results()[45].setName("Face_3_46")
Face_3.results()[46].setName("Face_3_47")
Face_3.results()[47].setName("Face_3_48")
Face_3.results()[48].setName("Face_3_49")
Face_3.results()[49].setName("Face_3_50")
Face_3.results()[50].setName("Face_3_51")
Face_3.results()[51].setName("Face_3_52")
Face_3.results()[52].setName("Face_3_53")
Face_3.results()[53].setName("Face_3_54")
Face_3.results()[54].setName("Face_3_55")
Face_3.results()[55].setName("Face_3_56")
Face_3.results()[56].setName("Face_3_57")
Face_3.results()[57].setName("Face_3_58")
Face_3.results()[58].setName("Face_3_59")
Face_3.results()[59].setName("Face_3_60")
Face_3.results()[60].setName("Face_3_61")
Face_3.results()[61].setName("Face_3_62")
Face_3.results()[62].setName("Face_3_63")
Face_3.results()[63].setName("Face_3_64")
Face_3.results()[64].setName("Face_3_65")
Face_3.results()[65].setName("Face_3_66")
Face_3.results()[66].setName("Face_3_67")
Face_3.results()[67].setName("Face_3_68")
Face_3.results()[68].setName("Face_3_69")
Face_3.results()[69].setName("Face_3_70")
Face_3.results()[70].setName("Face_3_71")

### Create Sketch
Sketch_4 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_9 = Sketch_4.addLine(0.677687, -0.0831312344139651, 0.677687, -0.07311876558603492)

### Create SketchLine
SketchLine_10 = Sketch_4.addLine(0.6801807655860349, -0.07062500000000001, 0.6801932344139652, -0.07062500000000663)

### Create SketchLine
SketchLine_11 = Sketch_4.addLine(0.682687, -0.07311876558604155, 0.682687, -0.0831312344139651)

### Create SketchLine
SketchLine_12 = Sketch_4.addLine(0.6801932344139652, -0.08562500000000001, 0.6801807655860349, -0.08562500000000001)

### Create SketchPoint
SketchPoint_7 = Sketch_4.addPoint(0.677687, -0.078125)
SketchPoint_7.setAuxiliary(True)
Sketch_4.setCoincident(SketchPoint_7.result(), SketchLine_9.result())
Sketch_4.setFixed(SketchPoint_7.result())

### Create SketchPoint
SketchPoint_8 = Sketch_4.addPoint(0.6801870000000001, -0.07062500000000001)
SketchPoint_8.setAuxiliary(True)
Sketch_4.setCoincident(SketchPoint_8.result(), SketchLine_10.result())
Sketch_4.setFixed(SketchPoint_8.result())

### Create SketchPoint
SketchPoint_9 = Sketch_4.addPoint(0.682687, -0.078125)
SketchPoint_9.setAuxiliary(True)
Sketch_4.setCoincident(SketchPoint_9.result(), SketchLine_11.result())
Sketch_4.setFixed(SketchPoint_9.result())

### Create SketchPoint
SketchPoint_10 = Sketch_4.addPoint(0.6801870000000001, -0.08562500000000001)
SketchPoint_10.setAuxiliary(True)
Sketch_4.setCoincident(SketchPoint_10.result(), SketchLine_12.result())
Sketch_4.setFixed(SketchPoint_10.result())

### Create SketchArc
SketchArc_5 = Sketch_4.addArc(0.6801807655860349, -0.0831312344139651, 0.677687, -0.0831312344139651, 0.6801807655860349, -0.08562500000000001, False)
Sketch_4.setCoincident(SketchArc_5.startPoint(), SketchLine_9.startPoint())
Sketch_4.setCoincident(SketchArc_5.endPoint(), SketchLine_12.endPoint())
Sketch_4.setTangent(SketchArc_5.results()[1], SketchLine_9.result())
Sketch_4.setTangent(SketchArc_5.results()[1], SketchLine_12.result())
Sketch_4.setRadius(SketchArc_5.results()[1], 0.002493765586034913)

### Create SketchArc
SketchArc_6 = Sketch_4.addArc(0.6801807655860349, -0.07311876558603492, 0.6801807655860349, -0.07062500000000001, 0.677687, -0.07311876558603492, False)
Sketch_4.setCoincident(SketchArc_6.startPoint(), SketchLine_10.startPoint())
Sketch_4.setCoincident(SketchArc_6.endPoint(), SketchLine_9.endPoint())
Sketch_4.setTangent(SketchArc_6.results()[1], SketchLine_9.result())
Sketch_4.setTangent(SketchArc_6.results()[1], SketchLine_10.result())
Sketch_4.setRadius(SketchArc_6.results()[1], 0.002493765586034913)

### Create SketchArc
SketchArc_7 = Sketch_4.addArc(0.6801932344139652, -0.07311876558604155, 0.682687, -0.07311876558604155, 0.6801932344139652, -0.07062500000000663, False)
Sketch_4.setCoincident(SketchArc_7.startPoint(), SketchLine_11.startPoint())
Sketch_4.setCoincident(SketchArc_7.endPoint(), SketchLine_10.endPoint())
Sketch_4.setTangent(SketchArc_7.results()[1], SketchLine_10.result())
Sketch_4.setTangent(SketchArc_7.results()[1], SketchLine_11.result())
Sketch_4.setRadius(SketchArc_7.results()[1], 0.002493765586034913)

### Create SketchArc
SketchArc_8 = Sketch_4.addArc(0.6801932344139652, -0.0831312344139651, 0.6801932344139652, -0.08562500000000001, 0.682687, -0.0831312344139651, False)
Sketch_4.setCoincident(SketchArc_8.startPoint(), SketchLine_12.startPoint())
Sketch_4.setCoincident(SketchArc_8.endPoint(), SketchLine_11.endPoint())
Sketch_4.setTangent(SketchArc_8.results()[1], SketchLine_12.result())
Sketch_4.setTangent(SketchArc_8.results()[1], SketchLine_11.result())
Sketch_4.setRadius(SketchArc_8.results()[1], 0.002493765586034913)

### Create SketchPoint
SketchPoint_11 = Sketch_4.addPoint(0.680187, -0.078125)

### Create SketchPoint
SketchPoint_12 = Sketch_4.addPoint(0.680187, -0.100125)

### Create SketchMultiTranslation
SketchMultiTranslation_4_objects = [SketchLine_9.result(), SketchLine_10.result(), SketchLine_11.result(), SketchLine_12.result(), SketchArc_5.results()[1], SketchArc_6.results()[1], SketchArc_7.results()[1], SketchArc_8.results()[1]]
SketchMultiTranslation_4 = Sketch_4.addTranslation(SketchMultiTranslation_4_objects, SketchPoint_11.coordinates(), SketchPoint_12.coordinates(), 64)
[SketchLine_13, SketchLine_14, SketchLine_15, SketchLine_16, SketchLine_17, SketchLine_18, SketchLine_19, SketchLine_20, SketchLine_21, SketchLine_22, SketchLine_23, SketchLine_24, SketchLine_25, SketchLine_26, SketchLine_27, SketchLine_28, SketchLine_29, SketchLine_30, SketchLine_31, SketchLine_32, SketchLine_33, SketchLine_34, SketchLine_35, SketchLine_36, SketchLine_37, SketchLine_38, SketchLine_39, SketchLine_40, SketchLine_41, SketchLine_42, SketchLine_43, SketchLine_44, SketchLine_45, SketchLine_46, SketchLine_47, SketchLine_48, SketchLine_49, SketchLine_50, SketchLine_51, SketchLine_52, SketchLine_53, SketchLine_54, SketchLine_55, SketchLine_56, SketchLine_57, SketchLine_58, SketchLine_59, SketchLine_60, SketchLine_61, SketchLine_62, SketchLine_63, SketchLine_64, SketchLine_65, SketchLine_66, SketchLine_67, SketchLine_68, SketchLine_69, SketchLine_70, SketchLine_71, SketchLine_72, SketchLine_73, SketchLine_74, SketchLine_75, SketchLine_76, SketchLine_77, SketchLine_78, SketchLine_79, SketchLine_80, SketchLine_81, SketchLine_82, SketchLine_83, SketchLine_84, SketchLine_85, SketchLine_86, SketchLine_87, SketchLine_88, SketchLine_89, SketchLine_90, SketchLine_91, SketchLine_92, SketchLine_93, SketchLine_94, SketchLine_95, SketchLine_96, SketchLine_97, SketchLine_98, SketchLine_99, SketchLine_100, SketchLine_101, SketchLine_102, SketchLine_103, SketchLine_104, SketchLine_105, SketchLine_106, SketchLine_107, SketchLine_108, SketchLine_109, SketchLine_110, SketchLine_111, SketchLine_112, SketchLine_113, SketchLine_114, SketchLine_115, SketchLine_116, SketchLine_117, SketchLine_118, SketchLine_119, SketchLine_120, SketchLine_121, SketchLine_122, SketchLine_123, SketchLine_124, SketchLine_125, SketchLine_126, SketchLine_127, SketchLine_128, SketchLine_129, SketchLine_130, SketchLine_131, SketchLine_132, SketchLine_133, SketchLine_134, SketchLine_135, SketchLine_136, SketchLine_137, SketchLine_138, SketchLine_139, SketchLine_140, SketchLine_141, SketchLine_142, SketchLine_143, SketchLine_144, SketchLine_145, SketchLine_146, SketchLine_147, SketchLine_148, SketchLine_149, SketchLine_150, SketchLine_151, SketchLine_152, SketchLine_153, SketchLine_154, SketchLine_155, SketchLine_156, SketchLine_157, SketchLine_158, SketchLine_159, SketchLine_160, SketchLine_161, SketchLine_162, SketchLine_163, SketchLine_164, SketchLine_165, SketchLine_166, SketchLine_167, SketchLine_168, SketchLine_169, SketchLine_170, SketchLine_171, SketchLine_172, SketchLine_173, SketchLine_174, SketchLine_175, SketchLine_176, SketchLine_177, SketchLine_178, SketchLine_179, SketchLine_180, SketchLine_181, SketchLine_182, SketchLine_183, SketchLine_184, SketchLine_185, SketchLine_186, SketchLine_187, SketchLine_188, SketchLine_189, SketchLine_190, SketchLine_191, SketchLine_192, SketchLine_193, SketchLine_194, SketchLine_195, SketchLine_196, SketchLine_197, SketchLine_198, SketchLine_199, SketchLine_200, SketchLine_201, SketchLine_202, SketchLine_203, SketchLine_204, SketchLine_205, SketchLine_206, SketchLine_207, SketchLine_208, SketchLine_209, SketchLine_210, SketchLine_211, SketchLine_212, SketchLine_213, SketchLine_214, SketchLine_215, SketchLine_216, SketchLine_217, SketchLine_218, SketchLine_219, SketchLine_220, SketchLine_221, SketchLine_222, SketchLine_223, SketchLine_224, SketchLine_225, SketchLine_226, SketchLine_227, SketchLine_228, SketchLine_229, SketchLine_230, SketchLine_231, SketchLine_232, SketchLine_233, SketchLine_234, SketchLine_235, SketchLine_236, SketchLine_237, SketchLine_238, SketchLine_239, SketchLine_240, SketchLine_241, SketchLine_242, SketchLine_243, SketchLine_244, SketchLine_245, SketchLine_246, SketchLine_247, SketchLine_248, SketchLine_249, SketchLine_250, SketchLine_251, SketchLine_252, SketchLine_253, SketchLine_254, SketchLine_255, SketchLine_256, SketchLine_257, SketchLine_258, SketchLine_259, SketchLine_260, SketchLine_261, SketchLine_262, SketchLine_263, SketchLine_264, SketchArc_9, SketchArc_10, SketchArc_11, SketchArc_12, SketchArc_13, SketchArc_14, SketchArc_15, SketchArc_16, SketchArc_17, SketchArc_18, SketchArc_19, SketchArc_20, SketchArc_21, SketchArc_22, SketchArc_23, SketchArc_24, SketchArc_25, SketchArc_26, SketchArc_27, SketchArc_28, SketchArc_29, SketchArc_30, SketchArc_31, SketchArc_32, SketchArc_33, SketchArc_34, SketchArc_35, SketchArc_36, SketchArc_37, SketchArc_38, SketchArc_39, SketchArc_40, SketchArc_41, SketchArc_42, SketchArc_43, SketchArc_44, SketchArc_45, SketchArc_46, SketchArc_47, SketchArc_48, SketchArc_49, SketchArc_50, SketchArc_51, SketchArc_52, SketchArc_53, SketchArc_54, SketchArc_55, SketchArc_56, SketchArc_57, SketchArc_58, SketchArc_59, SketchArc_60, SketchArc_61, SketchArc_62, SketchArc_63, SketchArc_64, SketchArc_65, SketchArc_66, SketchArc_67, SketchArc_68, SketchArc_69, SketchArc_70, SketchArc_71, SketchArc_72, SketchArc_73, SketchArc_74, SketchArc_75, SketchArc_76, SketchArc_77, SketchArc_78, SketchArc_79, SketchArc_80, SketchArc_81, SketchArc_82, SketchArc_83, SketchArc_84, SketchArc_85, SketchArc_86, SketchArc_87, SketchArc_88, SketchArc_89, SketchArc_90, SketchArc_91, SketchArc_92, SketchArc_93, SketchArc_94, SketchArc_95, SketchArc_96, SketchArc_97, SketchArc_98, SketchArc_99, SketchArc_100, SketchArc_101, SketchArc_102, SketchArc_103, SketchArc_104, SketchArc_105, SketchArc_106, SketchArc_107, SketchArc_108, SketchArc_109, SketchArc_110, SketchArc_111, SketchArc_112, SketchArc_113, SketchArc_114, SketchArc_115, SketchArc_116, SketchArc_117, SketchArc_118, SketchArc_119, SketchArc_120, SketchArc_121, SketchArc_122, SketchArc_123, SketchArc_124, SketchArc_125, SketchArc_126, SketchArc_127, SketchArc_128, SketchArc_129, SketchArc_130, SketchArc_131, SketchArc_132, SketchArc_133, SketchArc_134, SketchArc_135, SketchArc_136, SketchArc_137, SketchArc_138, SketchArc_139, SketchArc_140, SketchArc_141, SketchArc_142, SketchArc_143, SketchArc_144, SketchArc_145, SketchArc_146, SketchArc_147, SketchArc_148, SketchArc_149, SketchArc_150, SketchArc_151, SketchArc_152, SketchArc_153, SketchArc_154, SketchArc_155, SketchArc_156, SketchArc_157, SketchArc_158, SketchArc_159, SketchArc_160, SketchArc_161, SketchArc_162, SketchArc_163, SketchArc_164, SketchArc_165, SketchArc_166, SketchArc_167, SketchArc_168, SketchArc_169, SketchArc_170, SketchArc_171, SketchArc_172, SketchArc_173, SketchArc_174, SketchArc_175, SketchArc_176, SketchArc_177, SketchArc_178, SketchArc_179, SketchArc_180, SketchArc_181, SketchArc_182, SketchArc_183, SketchArc_184, SketchArc_185, SketchArc_186, SketchArc_187, SketchArc_188, SketchArc_189, SketchArc_190, SketchArc_191, SketchArc_192, SketchArc_193, SketchArc_194, SketchArc_195, SketchArc_196, SketchArc_197, SketchArc_198, SketchArc_199, SketchArc_200, SketchArc_201, SketchArc_202, SketchArc_203, SketchArc_204, SketchArc_205, SketchArc_206, SketchArc_207, SketchArc_208, SketchArc_209, SketchArc_210, SketchArc_211, SketchArc_212, SketchArc_213, SketchArc_214, SketchArc_215, SketchArc_216, SketchArc_217, SketchArc_218, SketchArc_219, SketchArc_220, SketchArc_221, SketchArc_222, SketchArc_223, SketchArc_224, SketchArc_225, SketchArc_226, SketchArc_227, SketchArc_228, SketchArc_229, SketchArc_230, SketchArc_231, SketchArc_232, SketchArc_233, SketchArc_234, SketchArc_235, SketchArc_236, SketchArc_237, SketchArc_238, SketchArc_239, SketchArc_240, SketchArc_241, SketchArc_242, SketchArc_243, SketchArc_244, SketchArc_245, SketchArc_246, SketchArc_247, SketchArc_248, SketchArc_249, SketchArc_250, SketchArc_251, SketchArc_252, SketchArc_253, SketchArc_254, SketchArc_255, SketchArc_256, SketchArc_257, SketchArc_258, SketchArc_259, SketchArc_260] = SketchMultiTranslation_4.translatedList()
model.do()
Sketch_4.setName("FieldCageGuard")
Sketch_4.result().setName("FieldCageGuard")

### Create Face
Face_4 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "FieldCageGuard")])
Face_4.setName("FieldCageGuard")
Face_4.result().setName("Face_4_1")
Face_4.results()[1].setName("Face_4_2")
Face_4.results()[2].setName("Face_4_3")
Face_4.results()[3].setName("Face_4_4")
Face_4.results()[4].setName("Face_4_5")
Face_4.results()[5].setName("Face_4_6")
Face_4.results()[6].setName("Face_4_7")
Face_4.results()[7].setName("Face_4_8")
Face_4.results()[8].setName("Face_4_9")
Face_4.results()[9].setName("Face_4_10")
Face_4.results()[10].setName("Face_4_11")
Face_4.results()[11].setName("Face_4_12")
Face_4.results()[12].setName("Face_4_13")
Face_4.results()[13].setName("Face_4_14")
Face_4.results()[14].setName("Face_4_15")
Face_4.results()[15].setName("Face_4_16")
Face_4.results()[16].setName("Face_4_17")
Face_4.results()[17].setName("Face_4_18")
Face_4.results()[18].setName("Face_4_19")
Face_4.results()[19].setName("Face_4_20")
Face_4.results()[20].setName("Face_4_21")
Face_4.results()[21].setName("Face_4_22")
Face_4.results()[22].setName("Face_4_23")
Face_4.results()[23].setName("Face_4_24")
Face_4.results()[24].setName("Face_4_25")
Face_4.results()[25].setName("Face_4_26")
Face_4.results()[26].setName("Face_4_27")
Face_4.results()[27].setName("Face_4_28")
Face_4.results()[28].setName("Face_4_29")
Face_4.results()[29].setName("Face_4_30")
Face_4.results()[30].setName("Face_4_31")
Face_4.results()[31].setName("Face_4_32")
Face_4.results()[32].setName("Face_4_33")
Face_4.results()[33].setName("Face_4_34")
Face_4.results()[34].setName("Face_4_35")
Face_4.results()[35].setName("Face_4_36")
Face_4.results()[36].setName("Face_4_37")
Face_4.results()[37].setName("Face_4_38")
Face_4.results()[38].setName("Face_4_39")
Face_4.results()[39].setName("Face_4_40")
Face_4.results()[40].setName("Face_4_41")
Face_4.results()[41].setName("Face_4_42")
Face_4.results()[42].setName("Face_4_43")
Face_4.results()[43].setName("Face_4_44")
Face_4.results()[44].setName("Face_4_45")
Face_4.results()[45].setName("Face_4_46")
Face_4.results()[46].setName("Face_4_47")
Face_4.results()[47].setName("Face_4_48")
Face_4.results()[48].setName("Face_4_49")
Face_4.results()[49].setName("Face_4_50")
Face_4.results()[50].setName("Face_4_51")
Face_4.results()[51].setName("Face_4_52")
Face_4.results()[52].setName("Face_4_53")
Face_4.results()[53].setName("Face_4_54")
Face_4.results()[54].setName("Face_4_55")
Face_4.results()[55].setName("Face_4_56")
Face_4.results()[56].setName("Face_4_57")
Face_4.results()[57].setName("Face_4_58")
Face_4.results()[58].setName("Face_4_59")
Face_4.results()[59].setName("Face_4_60")
Face_4.results()[60].setName("Face_4_61")
Face_4.results()[61].setName("Face_4_62")
Face_4.results()[62].setName("Face_4_63")
Face_4.results()[63].setName("Face_4_64")

### Create Sketch
Sketch_5 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_265 = Sketch_5.addLine(0, 0.239, 0.713, 0.239)

### Create SketchLine
SketchLine_266 = Sketch_5.addLine(0.713, 0.239, 0.713, -0.02500000000000002)

### Create SketchLine
SketchLine_267 = Sketch_5.addLine(0.713, -0.02500000000000002, 0.709, -0.02500000000000002)

### Create SketchLine
SketchLine_268 = Sketch_5.addLine(0.709, -0.02500000000000002, 0.709, 0.159)

### Create SketchLine
SketchLine_269 = Sketch_5.addLine(0.709, 0.159, 0.708, 0.159)

### Create SketchLine
SketchLine_270 = Sketch_5.addLine(0.708, 0.159, 0.708, 0.234)

### Create SketchLine
SketchLine_271 = Sketch_5.addLine(0.708, 0.234, 0, 0.234)

### Create SketchLine
SketchLine_272 = Sketch_5.addLine(0, 0.234, 0, 0.239)
Sketch_5.setCoincident(SketchLine_265.endPoint(), SketchLine_266.startPoint())
Sketch_5.setCoincident(SketchLine_266.endPoint(), SketchLine_267.startPoint())
Sketch_5.setCoincident(SketchLine_267.endPoint(), SketchLine_268.startPoint())
Sketch_5.setCoincident(SketchLine_268.endPoint(), SketchLine_269.startPoint())
Sketch_5.setCoincident(SketchLine_269.endPoint(), SketchLine_270.startPoint())
Sketch_5.setCoincident(SketchLine_270.endPoint(), SketchLine_271.startPoint())
Sketch_5.setCoincident(SketchLine_271.endPoint(), SketchLine_272.startPoint())
Sketch_5.setCoincident(SketchLine_272.endPoint(), SketchLine_265.startPoint())
model.do()
Sketch_5.setName("Bell")
Sketch_5.result().setName("Bell")

### Create Face
Face_5 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "Bell")])
Face_5.setName("Bell")
Face_5.result().setName("Face_5_1")

### Create Sketch
Sketch_6 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_273 = Sketch_6.addLine(0.667, -0.0074, 0.667, -0.001599999999999998)

### Create SketchLine
SketchLine_274 = Sketch_6.addLine(0.6686000000000001, 2.620744150175354e-18, 0.6964, -2.42887165707104e-18)

### Create SketchLine
SketchLine_275 = Sketch_6.addLine(0.6980000000000001, -0.001600000000000003, 0.6980000000000001, -0.01816666666666667)

### Create SketchLine
SketchLine_276 = Sketch_6.addLine(0.6964, -0.01976666666666667, 0.6788333333333334, -0.01976666666666667)

### Create SketchLine
SketchLine_277 = Sketch_6.addLine(0.6772333333333334, -0.01816666666666667, 0.6772333333333334, -0.009000000000000001)

### Create SketchLine
SketchLine_278 = Sketch_6.addLine(0.6772333333333334, -0.009000000000000001, 0.6686000000000001, -0.009000000000000001)
Sketch_6.setCoincident(SketchLine_277.endPoint(), SketchLine_278.startPoint())

### Create SketchPoint
SketchPoint_13 = Sketch_6.addPoint(0.667, -0.004500000000000001)
SketchPoint_13.setAuxiliary(True)
Sketch_6.setCoincident(SketchPoint_13.result(), SketchLine_273.result())
Sketch_6.setFixed(SketchPoint_13.result())

### Create SketchPoint
SketchPoint_14 = Sketch_6.addPoint(0.6825000000000001, 0)
SketchPoint_14.setAuxiliary(True)
Sketch_6.setCoincident(SketchPoint_14.result(), SketchLine_274.result())
Sketch_6.setFixed(SketchPoint_14.result())

### Create SketchArc
SketchArc_261 = Sketch_6.addArc(0.6686000000000001, -0.0074, 0.667, -0.0074, 0.6686000000000001, -0.009000000000000001, False)
Sketch_6.setCoincident(SketchArc_261.startPoint(), SketchLine_273.startPoint())
Sketch_6.setCoincident(SketchArc_261.endPoint(), SketchLine_278.endPoint())
Sketch_6.setTangent(SketchArc_261.results()[1], SketchLine_278.result())
Sketch_6.setTangent(SketchArc_261.results()[1], SketchLine_273.result())
Sketch_6.setRadius(SketchArc_261.results()[1], 0.0016)

### Create SketchArc
SketchArc_262 = Sketch_6.addArc(0.6686000000000001, -0.001599999999999997, 0.6686000000000001, 2.620744150175354e-18, 0.667, -0.001599999999999998, False)
Sketch_6.setCoincident(SketchArc_262.startPoint(), SketchLine_274.startPoint())
Sketch_6.setCoincident(SketchArc_262.endPoint(), SketchLine_273.endPoint())
Sketch_6.setTangent(SketchArc_262.results()[1], SketchLine_274.result())
Sketch_6.setTangent(SketchArc_262.results()[1], SketchLine_273.result())
Sketch_6.setRadius(SketchArc_262.results()[1], 0.0016)

### Create SketchArc
SketchArc_263 = Sketch_6.addArc(0.6964, -0.001600000000000003, 0.6980000000000001, -0.001600000000000003, 0.6964, -2.42887165707104e-18, False)
Sketch_6.setCoincident(SketchArc_263.startPoint(), SketchLine_275.startPoint())
Sketch_6.setCoincident(SketchArc_263.endPoint(), SketchLine_274.endPoint())
Sketch_6.setTangent(SketchArc_263.results()[1], SketchLine_275.result())
Sketch_6.setTangent(SketchArc_263.results()[1], SketchLine_274.result())
Sketch_6.setRadius(SketchArc_263.results()[1], 0.0016)

### Create SketchArc
SketchArc_264 = Sketch_6.addArc(0.6964, -0.01816666666666667, 0.6964, -0.01976666666666667, 0.6980000000000001, -0.01816666666666667, False)
Sketch_6.setCoincident(SketchArc_264.startPoint(), SketchLine_276.startPoint())
Sketch_6.setCoincident(SketchArc_264.endPoint(), SketchLine_275.endPoint())
Sketch_6.setTangent(SketchArc_264.results()[1], SketchLine_275.result())
Sketch_6.setTangent(SketchArc_264.results()[1], SketchLine_276.result())
Sketch_6.setRadius(SketchArc_264.results()[1], 0.0016)

### Create SketchArc
SketchArc_265 = Sketch_6.addArc(0.6788333333333334, -0.01816666666666667, 0.6772333333333334, -0.01816666666666667, 0.6788333333333334, -0.01976666666666667, False)
Sketch_6.setCoincident(SketchArc_265.startPoint(), SketchLine_277.startPoint())
Sketch_6.setCoincident(SketchArc_265.endPoint(), SketchLine_276.endPoint())
Sketch_6.setTangent(SketchArc_265.results()[1], SketchLine_276.result())
Sketch_6.setTangent(SketchArc_265.results()[1], SketchLine_277.result())
Sketch_6.setRadius(SketchArc_265.results()[1], 0.0016)

### Create SketchCircle
SketchCircle_72 = Sketch_6.addCircle(0.52625, 0, 5e-05)

### Create SketchCircle
SketchCircle_73 = Sketch_6.addCircle(0.24125, 0, 5e-05)

### Create SketchCircle
SketchCircle_74 = Sketch_6.addCircle(0.42125, 0, 5e-05)

### Create SketchCircle
SketchCircle_75 = Sketch_6.addCircle(0.50625, 0, 5e-05)

### Create SketchCircle
SketchCircle_76 = Sketch_6.addCircle(0.01625, 0, 5e-05)

### Create SketchCircle
SketchCircle_77 = Sketch_6.addCircle(0.14625, 0, 5e-05)

### Create SketchCircle
SketchCircle_78 = Sketch_6.addCircle(0.21625, 0, 5e-05)

### Create SketchCircle
SketchCircle_79 = Sketch_6.addCircle(0.38625, 0, 5e-05)

### Create SketchCircle
SketchCircle_80 = Sketch_6.addCircle(0.14125, 0, 5e-05)

### Create SketchCircle
SketchCircle_81 = Sketch_6.addCircle(0.36625, 0, 5e-05)

### Create SketchCircle
SketchCircle_82 = Sketch_6.addCircle(0.39125, 0, 5e-05)

### Create SketchCircle
SketchCircle_83 = Sketch_6.addCircle(0.08125, 0, 5e-05)

### Create SketchCircle
SketchCircle_84 = Sketch_6.addCircle(0.03625, 0, 5e-05)

### Create SketchCircle
SketchCircle_85 = Sketch_6.addCircle(0.08625000000000001, 0, 5e-05)

### Create SketchCircle
SketchCircle_86 = Sketch_6.addCircle(0.10625, 0, 5e-05)

### Create SketchCircle
SketchCircle_87 = Sketch_6.addCircle(0.07125000000000001, 0, 5e-05)

### Create SketchCircle
SketchCircle_88 = Sketch_6.addCircle(0.00625, 0, 5e-05)

### Create SketchCircle
SketchCircle_89 = Sketch_6.addCircle(0.63625, 0, 5e-05)

### Create SketchCircle
SketchCircle_90 = Sketch_6.addCircle(0.16625, 0, 5e-05)

### Create SketchCircle
SketchCircle_91 = Sketch_6.addCircle(0.66125, 0, 5e-05)

### Create SketchCircle
SketchCircle_92 = Sketch_6.addCircle(0.60125, 0, 5e-05)

### Create SketchCircle
SketchCircle_93 = Sketch_6.addCircle(0.21125, 0, 5e-05)

### Create SketchCircle
SketchCircle_94 = Sketch_6.addCircle(0.57125, 0, 5e-05)

### Create SketchCircle
SketchCircle_95 = Sketch_6.addCircle(0.26125, 0, 5e-05)

### Create SketchCircle
SketchCircle_96 = Sketch_6.addCircle(0.42625, 0, 5e-05)

### Create SketchCircle
SketchCircle_97 = Sketch_6.addCircle(0.02625, 0, 5e-05)

### Create SketchCircle
SketchCircle_98 = Sketch_6.addCircle(0.44125, 0, 5e-05)

### Create SketchCircle
SketchCircle_99 = Sketch_6.addCircle(0.19125, 0, 5e-05)

### Create SketchCircle
SketchCircle_100 = Sketch_6.addCircle(0.43625, 0, 5e-05)

### Create SketchCircle
SketchCircle_101 = Sketch_6.addCircle(0.25625, 0, 5e-05)

### Create SketchCircle
SketchCircle_102 = Sketch_6.addCircle(0.22625, 0, 5e-05)

### Create SketchCircle
SketchCircle_103 = Sketch_6.addCircle(0.07625, 0, 5e-05)

### Create SketchCircle
SketchCircle_104 = Sketch_6.addCircle(0.31125, 0, 5e-05)

### Create SketchCircle
SketchCircle_105 = Sketch_6.addCircle(0.51125, 0, 5e-05)

### Create SketchCircle
SketchCircle_106 = Sketch_6.addCircle(0.49125, 0, 5e-05)

### Create SketchCircle
SketchCircle_107 = Sketch_6.addCircle(0.30125, 0, 5e-05)

### Create SketchCircle
SketchCircle_108 = Sketch_6.addCircle(0.33625, 0, 5e-05)

### Create SketchCircle
SketchCircle_109 = Sketch_6.addCircle(0.45625, 0, 5e-05)

### Create SketchCircle
SketchCircle_110 = Sketch_6.addCircle(0.15125, 0, 5e-05)

### Create SketchCircle
SketchCircle_111 = Sketch_6.addCircle(0.36125, 0, 5e-05)

### Create SketchCircle
SketchCircle_112 = Sketch_6.addCircle(0.5912499999999999, 0, 5e-05)

### Create SketchCircle
SketchCircle_113 = Sketch_6.addCircle(0.04125, 0, 5e-05)

### Create SketchCircle
SketchCircle_114 = Sketch_6.addCircle(0.23125, 0, 5e-05)

### Create SketchCircle
SketchCircle_115 = Sketch_6.addCircle(0.32625, 0, 5e-05)

### Create SketchCircle
SketchCircle_116 = Sketch_6.addCircle(0.30625, 0, 5e-05)

### Create SketchCircle
SketchCircle_117 = Sketch_6.addCircle(0.37125, 0, 5e-05)

### Create SketchCircle
SketchCircle_118 = Sketch_6.addCircle(0.20625, 0, 5e-05)

### Create SketchCircle
SketchCircle_119 = Sketch_6.addCircle(0.25125, 0, 5e-05)

### Create SketchCircle
SketchCircle_120 = Sketch_6.addCircle(0.46125, 0, 5e-05)

### Create SketchCircle
SketchCircle_121 = Sketch_6.addCircle(0.29125, 0, 5e-05)

### Create SketchCircle
SketchCircle_122 = Sketch_6.addCircle(0.48125, 0, 5e-05)

### Create SketchCircle
SketchCircle_123 = Sketch_6.addCircle(0.46625, 0, 5e-05)

### Create SketchCircle
SketchCircle_124 = Sketch_6.addCircle(0.10125, 0, 5e-05)

### Create SketchCircle
SketchCircle_125 = Sketch_6.addCircle(0.31625, 0, 5e-05)

### Create SketchCircle
SketchCircle_126 = Sketch_6.addCircle(0.33125, 0, 5e-05)

### Create SketchCircle
SketchCircle_127 = Sketch_6.addCircle(0.5812499999999999, 0, 5e-05)

### Create SketchCircle
SketchCircle_128 = Sketch_6.addCircle(0.50125, 0, 5e-05)

### Create SketchCircle
SketchCircle_129 = Sketch_6.addCircle(0.13125, 0, 5e-05)

### Create SketchCircle
SketchCircle_130 = Sketch_6.addCircle(0.34125, 0, 5e-05)

### Create SketchCircle
SketchCircle_131 = Sketch_6.addCircle(0.64625, 0, 5e-05)

### Create SketchCircle
SketchCircle_132 = Sketch_6.addCircle(0.02125, 0, 5e-05)

### Create SketchCircle
SketchCircle_133 = Sketch_6.addCircle(0.17625, 0, 5e-05)

### Create SketchCircle
SketchCircle_134 = Sketch_6.addCircle(0.03125, 0, 5e-05)

### Create SketchCircle
SketchCircle_135 = Sketch_6.addCircle(0.04625, 0, 5e-05)

### Create SketchCircle
SketchCircle_136 = Sketch_6.addCircle(0.05625, 0, 5e-05)

### Create SketchCircle
SketchCircle_137 = Sketch_6.addCircle(0.00125, 0, 5e-05)

### Create SketchCircle
SketchCircle_138 = Sketch_6.addCircle(0.54625, 0, 5e-05)

### Create SketchCircle
SketchCircle_139 = Sketch_6.addCircle(0.5862499999999999, 0, 5e-05)

### Create SketchCircle
SketchCircle_140 = Sketch_6.addCircle(0.62125, 0, 5e-05)

### Create SketchCircle
SketchCircle_141 = Sketch_6.addCircle(0.40125, 0, 5e-05)

### Create SketchCircle
SketchCircle_142 = Sketch_6.addCircle(0.22125, 0, 5e-05)

### Create SketchCircle
SketchCircle_143 = Sketch_6.addCircle(0.23625, 0, 5e-05)

### Create SketchCircle
SketchCircle_144 = Sketch_6.addCircle(0.51625, 0, 5e-05)

### Create SketchCircle
SketchCircle_145 = Sketch_6.addCircle(0.65125, 0, 5e-05)

### Create SketchCircle
SketchCircle_146 = Sketch_6.addCircle(0.61125, 0, 5e-05)

### Create SketchCircle
SketchCircle_147 = Sketch_6.addCircle(0.11625, 0, 5e-05)

### Create SketchCircle
SketchCircle_148 = Sketch_6.addCircle(0.65625, 0, 5e-05)

### Create SketchCircle
SketchCircle_149 = Sketch_6.addCircle(0.15625, 0, 5e-05)

### Create SketchCircle
SketchCircle_150 = Sketch_6.addCircle(0.55125, 0, 5e-05)

### Create SketchCircle
SketchCircle_151 = Sketch_6.addCircle(0.24625, 0, 5e-05)

### Create SketchCircle
SketchCircle_152 = Sketch_6.addCircle(0.48625, 0, 5e-05)

### Create SketchCircle
SketchCircle_153 = Sketch_6.addCircle(0.47625, 0, 5e-05)

### Create SketchCircle
SketchCircle_154 = Sketch_6.addCircle(0.16125, 0, 5e-05)

### Create SketchCircle
SketchCircle_155 = Sketch_6.addCircle(0.05125, 0, 5e-05)

### Create SketchCircle
SketchCircle_156 = Sketch_6.addCircle(0.38125, 0, 5e-05)

### Create SketchCircle
SketchCircle_157 = Sketch_6.addCircle(0.37625, 0, 5e-05)

### Create SketchCircle
SketchCircle_158 = Sketch_6.addCircle(0.18625, 0, 5e-05)

### Create SketchCircle
SketchCircle_159 = Sketch_6.addCircle(0.19625, 0, 5e-05)

### Create SketchCircle
SketchCircle_160 = Sketch_6.addCircle(0.28125, 0, 5e-05)

### Create SketchCircle
SketchCircle_161 = Sketch_6.addCircle(0.53125, 0, 5e-05)

### Create SketchCircle
SketchCircle_162 = Sketch_6.addCircle(0.34625, 0, 5e-05)

### Create SketchCircle
SketchCircle_163 = Sketch_6.addCircle(0.01125, 0, 5e-05)

### Create SketchCircle
SketchCircle_164 = Sketch_6.addCircle(0.09125, 0, 5e-05)

### Create SketchCircle
SketchCircle_165 = Sketch_6.addCircle(0.17125, 0, 5e-05)

### Create SketchCircle
SketchCircle_166 = Sketch_6.addCircle(0.56125, 0, 5e-05)

### Create SketchCircle
SketchCircle_167 = Sketch_6.addCircle(0.53625, 0, 5e-05)

### Create SketchCircle
SketchCircle_168 = Sketch_6.addCircle(0.27625, 0, 5e-05)

### Create SketchCircle
SketchCircle_169 = Sketch_6.addCircle(0.5962499999999999, 0, 5e-05)

### Create SketchCircle
SketchCircle_170 = Sketch_6.addCircle(0.40625, 0, 5e-05)

### Create SketchCircle
SketchCircle_171 = Sketch_6.addCircle(0.52125, 0, 5e-05)

### Create SketchCircle
SketchCircle_172 = Sketch_6.addCircle(0.49625, 0, 5e-05)

### Create SketchCircle
SketchCircle_173 = Sketch_6.addCircle(0.13625, 0, 5e-05)

### Create SketchCircle
SketchCircle_174 = Sketch_6.addCircle(0.26625, 0, 5e-05)

### Create SketchCircle
SketchCircle_175 = Sketch_6.addCircle(0.39625, 0, 5e-05)

### Create SketchCircle
SketchCircle_176 = Sketch_6.addCircle(0.06125, 0, 5e-05)

### Create SketchCircle
SketchCircle_177 = Sketch_6.addCircle(0.56625, 0, 5e-05)

### Create SketchCircle
SketchCircle_178 = Sketch_6.addCircle(0.41125, 0, 5e-05)

### Create SketchCircle
SketchCircle_179 = Sketch_6.addCircle(0.20125, 0, 5e-05)

### Create SketchCircle
SketchCircle_180 = Sketch_6.addCircle(0.62625, 0, 5e-05)

### Create SketchCircle
SketchCircle_181 = Sketch_6.addCircle(0.57625, 0, 5e-05)

### Create SketchCircle
SketchCircle_182 = Sketch_6.addCircle(0.55625, 0, 5e-05)

### Create SketchCircle
SketchCircle_183 = Sketch_6.addCircle(0.47125, 0, 5e-05)

### Create SketchCircle
SketchCircle_184 = Sketch_6.addCircle(0.63125, 0, 5e-05)

### Create SketchCircle
SketchCircle_185 = Sketch_6.addCircle(0.32125, 0, 5e-05)

### Create SketchCircle
SketchCircle_186 = Sketch_6.addCircle(0.09625, 0, 5e-05)

### Create SketchCircle
SketchCircle_187 = Sketch_6.addCircle(0.60625, 0, 5e-05)

### Create SketchCircle
SketchCircle_188 = Sketch_6.addCircle(0.64125, 0, 5e-05)

### Create SketchCircle
SketchCircle_189 = Sketch_6.addCircle(0.45125, 0, 5e-05)

### Create SketchCircle
SketchCircle_190 = Sketch_6.addCircle(0.28625, 0, 5e-05)

### Create SketchCircle
SketchCircle_191 = Sketch_6.addCircle(0.27125, 0, 5e-05)

### Create SketchCircle
SketchCircle_192 = Sketch_6.addCircle(0.41625, 0, 5e-05)

### Create SketchCircle
SketchCircle_193 = Sketch_6.addCircle(0.18125, 0, 5e-05)

### Create SketchCircle
SketchCircle_194 = Sketch_6.addCircle(0.11125, 0, 5e-05)

### Create SketchCircle
SketchCircle_195 = Sketch_6.addCircle(0.43125, 0, 5e-05)

### Create SketchCircle
SketchCircle_196 = Sketch_6.addCircle(0.06625, 0, 5e-05)

### Create SketchCircle
SketchCircle_197 = Sketch_6.addCircle(0.35125, 0, 5e-05)

### Create SketchCircle
SketchCircle_198 = Sketch_6.addCircle(0.12625, 0, 5e-05)

### Create SketchCircle
SketchCircle_199 = Sketch_6.addCircle(0.54125, 0, 5e-05)

### Create SketchCircle
SketchCircle_200 = Sketch_6.addCircle(0.12125, 0, 5e-05)

### Create SketchCircle
SketchCircle_201 = Sketch_6.addCircle(0.35625, 0, 5e-05)

### Create SketchCircle
SketchCircle_202 = Sketch_6.addCircle(0.29625, 0, 5e-05)

### Create SketchCircle
SketchCircle_203 = Sketch_6.addCircle(0.61625, 0, 5e-05)

### Create SketchCircle
SketchCircle_204 = Sketch_6.addCircle(0.44625, 0, 5e-05)
model.do()
Sketch_6.setName("Gate")
Sketch_6.result().setName("Gate")

### Create Face
Face_6 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "Gate")])
Face_6.setName("Gate")
Face_6.result().setName("Face_6_1")
Face_6.results()[1].setName("Face_6_2")
Face_6.results()[2].setName("Face_6_3")
Face_6.results()[3].setName("Face_6_4")
Face_6.results()[4].setName("Face_6_5")
Face_6.results()[5].setName("Face_6_6")
Face_6.results()[6].setName("Face_6_7")
Face_6.results()[7].setName("Face_6_8")
Face_6.results()[8].setName("Face_6_9")
Face_6.results()[9].setName("Face_6_10")
Face_6.results()[10].setName("Face_6_11")
Face_6.results()[11].setName("Face_6_12")
Face_6.results()[12].setName("Face_6_13")
Face_6.results()[13].setName("Face_6_14")
Face_6.results()[14].setName("Face_6_15")
Face_6.results()[15].setName("Face_6_16")
Face_6.results()[16].setName("Face_6_17")
Face_6.results()[17].setName("Face_6_18")
Face_6.results()[18].setName("Face_6_19")
Face_6.results()[19].setName("Face_6_20")
Face_6.results()[20].setName("Face_6_21")
Face_6.results()[21].setName("Face_6_22")
Face_6.results()[22].setName("Face_6_23")
Face_6.results()[23].setName("Face_6_24")
Face_6.results()[24].setName("Face_6_25")
Face_6.results()[25].setName("Face_6_26")
Face_6.results()[26].setName("Face_6_27")
Face_6.results()[27].setName("Face_6_28")
Face_6.results()[28].setName("Face_6_29")
Face_6.results()[29].setName("Face_6_30")
Face_6.results()[30].setName("Face_6_31")
Face_6.results()[31].setName("Face_6_32")
Face_6.results()[32].setName("Face_6_33")
Face_6.results()[33].setName("Face_6_34")
Face_6.results()[34].setName("Face_6_35")
Face_6.results()[35].setName("Face_6_36")
Face_6.results()[36].setName("Face_6_37")
Face_6.results()[37].setName("Face_6_38")
Face_6.results()[38].setName("Face_6_39")
Face_6.results()[39].setName("Face_6_40")
Face_6.results()[40].setName("Face_6_41")
Face_6.results()[41].setName("Face_6_42")
Face_6.results()[42].setName("Face_6_43")
Face_6.results()[43].setName("Face_6_44")
Face_6.results()[44].setName("Face_6_45")
Face_6.results()[45].setName("Face_6_46")
Face_6.results()[46].setName("Face_6_47")
Face_6.results()[47].setName("Face_6_48")
Face_6.results()[48].setName("Face_6_49")
Face_6.results()[49].setName("Face_6_50")
Face_6.results()[50].setName("Face_6_51")
Face_6.results()[51].setName("Face_6_52")
Face_6.results()[52].setName("Face_6_53")
Face_6.results()[53].setName("Face_6_54")
Face_6.results()[54].setName("Face_6_55")
Face_6.results()[55].setName("Face_6_56")
Face_6.results()[56].setName("Face_6_57")
Face_6.results()[57].setName("Face_6_58")
Face_6.results()[58].setName("Face_6_59")
Face_6.results()[59].setName("Face_6_60")
Face_6.results()[60].setName("Face_6_61")
Face_6.results()[61].setName("Face_6_62")
Face_6.results()[62].setName("Face_6_63")
Face_6.results()[63].setName("Face_6_64")
Face_6.results()[64].setName("Face_6_65")
Face_6.results()[65].setName("Face_6_66")
Face_6.results()[66].setName("Face_6_67")
Face_6.results()[67].setName("Face_6_68")
Face_6.results()[68].setName("Face_6_69")
Face_6.results()[69].setName("Face_6_70")
Face_6.results()[70].setName("Face_6_71")
Face_6.results()[71].setName("Face_6_72")
Face_6.results()[72].setName("Face_6_73")
Face_6.results()[73].setName("Face_6_74")
Face_6.results()[74].setName("Face_6_75")
Face_6.results()[75].setName("Face_6_76")
Face_6.results()[76].setName("Face_6_77")
Face_6.results()[77].setName("Face_6_78")
Face_6.results()[78].setName("Face_6_79")
Face_6.results()[79].setName("Face_6_80")
Face_6.results()[80].setName("Face_6_81")
Face_6.results()[81].setName("Face_6_82")
Face_6.results()[82].setName("Face_6_83")
Face_6.results()[83].setName("Face_6_84")
Face_6.results()[84].setName("Face_6_85")
Face_6.results()[85].setName("Face_6_86")
Face_6.results()[86].setName("Face_6_87")
Face_6.results()[87].setName("Face_6_88")
Face_6.results()[88].setName("Face_6_89")
Face_6.results()[89].setName("Face_6_90")
Face_6.results()[90].setName("Face_6_91")
Face_6.results()[91].setName("Face_6_92")
Face_6.results()[92].setName("Face_6_93")
Face_6.results()[93].setName("Face_6_94")
Face_6.results()[94].setName("Face_6_95")
Face_6.results()[95].setName("Face_6_96")
Face_6.results()[96].setName("Face_6_97")
Face_6.results()[97].setName("Face_6_98")
Face_6.results()[98].setName("Face_6_99")
Face_6.results()[99].setName("Face_6_100")
Face_6.results()[100].setName("Face_6_101")
Face_6.results()[101].setName("Face_6_102")
Face_6.results()[102].setName("Face_6_103")
Face_6.results()[103].setName("Face_6_104")
Face_6.results()[104].setName("Face_6_105")
Face_6.results()[105].setName("Face_6_106")
Face_6.results()[106].setName("Face_6_107")
Face_6.results()[107].setName("Face_6_108")
Face_6.results()[108].setName("Face_6_109")
Face_6.results()[109].setName("Face_6_110")
Face_6.results()[110].setName("Face_6_111")
Face_6.results()[111].setName("Face_6_112")
Face_6.results()[112].setName("Face_6_113")
Face_6.results()[113].setName("Face_6_114")
Face_6.results()[114].setName("Face_6_115")
Face_6.results()[115].setName("Face_6_116")
Face_6.results()[116].setName("Face_6_117")
Face_6.results()[117].setName("Face_6_118")
Face_6.results()[118].setName("Face_6_119")
Face_6.results()[119].setName("Face_6_120")
Face_6.results()[120].setName("Face_6_121")
Face_6.results()[121].setName("Face_6_122")
Face_6.results()[122].setName("Face_6_123")
Face_6.results()[123].setName("Face_6_124")
Face_6.results()[124].setName("Face_6_125")
Face_6.results()[125].setName("Face_6_126")
Face_6.results()[126].setName("Face_6_127")
Face_6.results()[127].setName("Face_6_128")
Face_6.results()[128].setName("Face_6_129")
Face_6.results()[129].setName("Face_6_130")
Face_6.results()[130].setName("Face_6_131")
Face_6.results()[131].setName("Face_6_132")
Face_6.results()[132].setName("Face_6_133")
Face_6.results()[133].setName("Face_6_134")

### Create Sketch
Sketch_7 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_279 = Sketch_7.addLine(0.667, 0.009600000000000001, 0.667, 0.0304)

### Create SketchLine
SketchLine_280 = Sketch_7.addLine(0.6686000000000001, 0.032, 0.6964, 0.032)

### Create SketchLine
SketchLine_281 = Sketch_7.addLine(0.6980000000000001, 0.0304, 0.6980000000000001, 0.009599999999999997)

### Create SketchLine
SketchLine_282 = Sketch_7.addLine(0.6964, 0.007999999999999998, 0.6686000000000001, 0.008)

### Create SketchPoint
SketchPoint_15 = Sketch_7.addPoint(0.667, 0.02)
SketchPoint_15.setAuxiliary(True)
Sketch_7.setCoincident(SketchPoint_15.result(), SketchLine_279.result())
Sketch_7.setFixed(SketchPoint_15.result())

### Create SketchPoint
SketchPoint_16 = Sketch_7.addPoint(0.6825000000000001, 0.032)
SketchPoint_16.setAuxiliary(True)
Sketch_7.setCoincident(SketchPoint_16.result(), SketchLine_280.result())
Sketch_7.setFixed(SketchPoint_16.result())

### Create SketchPoint
SketchPoint_17 = Sketch_7.addPoint(0.6980000000000001, 0.02)
SketchPoint_17.setAuxiliary(True)
Sketch_7.setCoincident(SketchPoint_17.result(), SketchLine_281.result())
Sketch_7.setFixed(SketchPoint_17.result())

### Create SketchPoint
SketchPoint_18 = Sketch_7.addPoint(0.6825000000000001, 0.008)
SketchPoint_18.setAuxiliary(True)
Sketch_7.setCoincident(SketchPoint_18.result(), SketchLine_282.result())
Sketch_7.setFixed(SketchPoint_18.result())

### Create SketchArc
SketchArc_266 = Sketch_7.addArc(0.6686000000000001, 0.009600000000000001, 0.667, 0.009600000000000001, 0.6686000000000001, 0.008, False)
Sketch_7.setCoincident(SketchArc_266.startPoint(), SketchLine_279.startPoint())
Sketch_7.setCoincident(SketchArc_266.endPoint(), SketchLine_282.endPoint())
Sketch_7.setTangent(SketchArc_266.results()[1], SketchLine_279.result())
Sketch_7.setTangent(SketchArc_266.results()[1], SketchLine_282.result())
Sketch_7.setRadius(SketchArc_266.results()[1], 0.0016)

### Create SketchArc
SketchArc_267 = Sketch_7.addArc(0.6686000000000001, 0.0304, 0.6686000000000001, 0.032, 0.667, 0.0304, False)
Sketch_7.setCoincident(SketchArc_267.startPoint(), SketchLine_280.startPoint())
Sketch_7.setCoincident(SketchArc_267.endPoint(), SketchLine_279.endPoint())
Sketch_7.setTangent(SketchArc_267.results()[1], SketchLine_280.result())
Sketch_7.setTangent(SketchArc_267.results()[1], SketchLine_279.result())
Sketch_7.setRadius(SketchArc_267.results()[1], 0.0016)

### Create SketchArc
SketchArc_268 = Sketch_7.addArc(0.6964, 0.0304, 0.6980000000000001, 0.0304, 0.6964, 0.032, False)
Sketch_7.setCoincident(SketchArc_268.startPoint(), SketchLine_281.startPoint())
Sketch_7.setCoincident(SketchArc_268.endPoint(), SketchLine_280.endPoint())
Sketch_7.setTangent(SketchArc_268.results()[1], SketchLine_280.result())
Sketch_7.setTangent(SketchArc_268.results()[1], SketchLine_281.result())
Sketch_7.setRadius(SketchArc_268.results()[1], 0.0016)

### Create SketchArc
SketchArc_269 = Sketch_7.addArc(0.6964, 0.009599999999999997, 0.6964, 0.007999999999999998, 0.6980000000000001, 0.009599999999999997, False)
Sketch_7.setCoincident(SketchArc_269.startPoint(), SketchLine_282.startPoint())
Sketch_7.setCoincident(SketchArc_269.endPoint(), SketchLine_281.endPoint())
Sketch_7.setTangent(SketchArc_269.results()[1], SketchLine_281.result())
Sketch_7.setTangent(SketchArc_269.results()[1], SketchLine_282.result())
Sketch_7.setRadius(SketchArc_269.results()[1], 0.0016)

### Create SketchCircle
SketchCircle_205 = Sketch_7.addCircle(0.47625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_206 = Sketch_7.addCircle(0.33125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_207 = Sketch_7.addCircle(0.49625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_208 = Sketch_7.addCircle(0.51125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_209 = Sketch_7.addCircle(0.06125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_210 = Sketch_7.addCircle(0.56125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_211 = Sketch_7.addCircle(0.16625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_212 = Sketch_7.addCircle(0.53625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_213 = Sketch_7.addCircle(0.01625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_214 = Sketch_7.addCircle(0.05625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_215 = Sketch_7.addCircle(0.5862499999999999, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_216 = Sketch_7.addCircle(0.51625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_217 = Sketch_7.addCircle(0.31125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_218 = Sketch_7.addCircle(0.36125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_219 = Sketch_7.addCircle(0.55625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_220 = Sketch_7.addCircle(0.04625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_221 = Sketch_7.addCircle(0.07125000000000001, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_222 = Sketch_7.addCircle(0.09125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_223 = Sketch_7.addCircle(0.12125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_224 = Sketch_7.addCircle(0.46125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_225 = Sketch_7.addCircle(0.32125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_226 = Sketch_7.addCircle(0.35125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_227 = Sketch_7.addCircle(0.41125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_228 = Sketch_7.addCircle(0.55125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_229 = Sketch_7.addCircle(0.54125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_230 = Sketch_7.addCircle(0.32625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_231 = Sketch_7.addCircle(0.24125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_232 = Sketch_7.addCircle(0.46625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_233 = Sketch_7.addCircle(0.30625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_234 = Sketch_7.addCircle(0.09625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_235 = Sketch_7.addCircle(0.11125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_236 = Sketch_7.addCircle(0.47125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_237 = Sketch_7.addCircle(0.49125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_238 = Sketch_7.addCircle(0.08625000000000001, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_239 = Sketch_7.addCircle(0.40125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_240 = Sketch_7.addCircle(0.39625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_241 = Sketch_7.addCircle(0.37125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_242 = Sketch_7.addCircle(0.03125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_243 = Sketch_7.addCircle(0.34125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_244 = Sketch_7.addCircle(0.18125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_245 = Sketch_7.addCircle(0.13625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_246 = Sketch_7.addCircle(0.22125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_247 = Sketch_7.addCircle(0.25125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_248 = Sketch_7.addCircle(0.33625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_249 = Sketch_7.addCircle(0.26625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_250 = Sketch_7.addCircle(0.14125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_251 = Sketch_7.addCircle(0.28625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_252 = Sketch_7.addCircle(0.26125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_253 = Sketch_7.addCircle(0.29625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_254 = Sketch_7.addCircle(0.04125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_255 = Sketch_7.addCircle(0.61625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_256 = Sketch_7.addCircle(0.17625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_257 = Sketch_7.addCircle(0.29125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_258 = Sketch_7.addCircle(0.34625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_259 = Sketch_7.addCircle(0.63625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_260 = Sketch_7.addCircle(0.56625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_261 = Sketch_7.addCircle(0.62625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_262 = Sketch_7.addCircle(0.31625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_263 = Sketch_7.addCircle(0.02125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_264 = Sketch_7.addCircle(0.03625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_265 = Sketch_7.addCircle(0.50125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_266 = Sketch_7.addCircle(0.19125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_267 = Sketch_7.addCircle(0.17125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_268 = Sketch_7.addCircle(0.48625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_269 = Sketch_7.addCircle(0.23625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_270 = Sketch_7.addCircle(0.60125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_271 = Sketch_7.addCircle(0.54625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_272 = Sketch_7.addCircle(0.01125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_273 = Sketch_7.addCircle(0.61125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_274 = Sketch_7.addCircle(0.13125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_275 = Sketch_7.addCircle(0.38125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_276 = Sketch_7.addCircle(0.15125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_277 = Sketch_7.addCircle(0.39125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_278 = Sketch_7.addCircle(0.40625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_279 = Sketch_7.addCircle(0.21625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_280 = Sketch_7.addCircle(0.28125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_281 = Sketch_7.addCircle(0.21125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_282 = Sketch_7.addCircle(0.16125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_283 = Sketch_7.addCircle(0.15625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_284 = Sketch_7.addCircle(0.14625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_285 = Sketch_7.addCircle(0.12625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_286 = Sketch_7.addCircle(0.06625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_287 = Sketch_7.addCircle(0.50625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_288 = Sketch_7.addCircle(0.37625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_289 = Sketch_7.addCircle(0.10125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_290 = Sketch_7.addCircle(0.57125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_291 = Sketch_7.addCircle(0.60625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_292 = Sketch_7.addCircle(0.11625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_293 = Sketch_7.addCircle(0.5962499999999999, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_294 = Sketch_7.addCircle(0.20125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_295 = Sketch_7.addCircle(0.42125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_296 = Sketch_7.addCircle(0.10625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_297 = Sketch_7.addCircle(0.27625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_298 = Sketch_7.addCircle(0.63125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_299 = Sketch_7.addCircle(0.35625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_300 = Sketch_7.addCircle(0.00625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_301 = Sketch_7.addCircle(0.53125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_302 = Sketch_7.addCircle(0.30125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_303 = Sketch_7.addCircle(0.23125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_304 = Sketch_7.addCircle(0.05125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_305 = Sketch_7.addCircle(0.64125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_306 = Sketch_7.addCircle(0.43125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_307 = Sketch_7.addCircle(0.45125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_308 = Sketch_7.addCircle(0.52125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_309 = Sketch_7.addCircle(0.52625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_310 = Sketch_7.addCircle(0.5812499999999999, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_311 = Sketch_7.addCircle(0.22625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_312 = Sketch_7.addCircle(0.41625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_313 = Sketch_7.addCircle(0.42625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_314 = Sketch_7.addCircle(0.43625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_315 = Sketch_7.addCircle(0.44125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_316 = Sketch_7.addCircle(0.08125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_317 = Sketch_7.addCircle(0.38625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_318 = Sketch_7.addCircle(0.62125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_319 = Sketch_7.addCircle(0.57625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_320 = Sketch_7.addCircle(0.20625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_321 = Sketch_7.addCircle(0.24625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_322 = Sketch_7.addCircle(0.66125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_323 = Sketch_7.addCircle(0.45625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_324 = Sketch_7.addCircle(0.48125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_325 = Sketch_7.addCircle(0.44625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_326 = Sketch_7.addCircle(0.36625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_327 = Sketch_7.addCircle(0.65125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_328 = Sketch_7.addCircle(0.64625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_329 = Sketch_7.addCircle(0.5912499999999999, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_330 = Sketch_7.addCircle(0.25625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_331 = Sketch_7.addCircle(0.07625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_332 = Sketch_7.addCircle(0.27125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_333 = Sketch_7.addCircle(0.19625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_334 = Sketch_7.addCircle(0.18625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_335 = Sketch_7.addCircle(0.00125, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_336 = Sketch_7.addCircle(0.65625, 0.008, 0.0001080000000000005)

### Create SketchCircle
SketchCircle_337 = Sketch_7.addCircle(0.02625, 0.008, 0.0001080000000000005)
model.do()
Sketch_7.setName("Anode")
Sketch_7.result().setName("Anode")

### Create Face
Face_7 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "Anode")])
Face_7.setName("Anode")
Face_7.result().setName("Face_7_1")
Face_7.results()[1].setName("Face_7_2")
Face_7.results()[2].setName("Face_7_3")
Face_7.results()[3].setName("Face_7_4")
Face_7.results()[4].setName("Face_7_5")
Face_7.results()[5].setName("Face_7_6")
Face_7.results()[6].setName("Face_7_7")
Face_7.results()[7].setName("Face_7_8")
Face_7.results()[8].setName("Face_7_9")
Face_7.results()[9].setName("Face_7_10")
Face_7.results()[10].setName("Face_7_11")
Face_7.results()[11].setName("Face_7_12")
Face_7.results()[12].setName("Face_7_13")
Face_7.results()[13].setName("Face_7_14")
Face_7.results()[14].setName("Face_7_15")
Face_7.results()[15].setName("Face_7_16")
Face_7.results()[16].setName("Face_7_17")
Face_7.results()[17].setName("Face_7_18")
Face_7.results()[18].setName("Face_7_19")
Face_7.results()[19].setName("Face_7_20")
Face_7.results()[20].setName("Face_7_21")
Face_7.results()[21].setName("Face_7_22")
Face_7.results()[22].setName("Face_7_23")
Face_7.results()[23].setName("Face_7_24")
Face_7.results()[24].setName("Face_7_25")
Face_7.results()[25].setName("Face_7_26")
Face_7.results()[26].setName("Face_7_27")
Face_7.results()[27].setName("Face_7_28")
Face_7.results()[28].setName("Face_7_29")
Face_7.results()[29].setName("Face_7_30")
Face_7.results()[30].setName("Face_7_31")
Face_7.results()[31].setName("Face_7_32")
Face_7.results()[32].setName("Face_7_33")
Face_7.results()[33].setName("Face_7_34")
Face_7.results()[34].setName("Face_7_35")
Face_7.results()[35].setName("Face_7_36")
Face_7.results()[36].setName("Face_7_37")
Face_7.results()[37].setName("Face_7_38")
Face_7.results()[38].setName("Face_7_39")
Face_7.results()[39].setName("Face_7_40")
Face_7.results()[40].setName("Face_7_41")
Face_7.results()[41].setName("Face_7_42")
Face_7.results()[42].setName("Face_7_43")
Face_7.results()[43].setName("Face_7_44")
Face_7.results()[44].setName("Face_7_45")
Face_7.results()[45].setName("Face_7_46")
Face_7.results()[46].setName("Face_7_47")
Face_7.results()[47].setName("Face_7_48")
Face_7.results()[48].setName("Face_7_49")
Face_7.results()[49].setName("Face_7_50")
Face_7.results()[50].setName("Face_7_51")
Face_7.results()[51].setName("Face_7_52")
Face_7.results()[52].setName("Face_7_53")
Face_7.results()[53].setName("Face_7_54")
Face_7.results()[54].setName("Face_7_55")
Face_7.results()[55].setName("Face_7_56")
Face_7.results()[56].setName("Face_7_57")
Face_7.results()[57].setName("Face_7_58")
Face_7.results()[58].setName("Face_7_59")
Face_7.results()[59].setName("Face_7_60")
Face_7.results()[60].setName("Face_7_61")
Face_7.results()[61].setName("Face_7_62")
Face_7.results()[62].setName("Face_7_63")
Face_7.results()[63].setName("Face_7_64")
Face_7.results()[64].setName("Face_7_65")
Face_7.results()[65].setName("Face_7_66")
Face_7.results()[66].setName("Face_7_67")
Face_7.results()[67].setName("Face_7_68")
Face_7.results()[68].setName("Face_7_69")
Face_7.results()[69].setName("Face_7_70")
Face_7.results()[70].setName("Face_7_71")
Face_7.results()[71].setName("Face_7_72")
Face_7.results()[72].setName("Face_7_73")
Face_7.results()[73].setName("Face_7_74")
Face_7.results()[74].setName("Face_7_75")
Face_7.results()[75].setName("Face_7_76")
Face_7.results()[76].setName("Face_7_77")
Face_7.results()[77].setName("Face_7_78")
Face_7.results()[78].setName("Face_7_79")
Face_7.results()[79].setName("Face_7_80")
Face_7.results()[80].setName("Face_7_81")
Face_7.results()[81].setName("Face_7_82")
Face_7.results()[82].setName("Face_7_83")
Face_7.results()[83].setName("Face_7_84")
Face_7.results()[84].setName("Face_7_85")
Face_7.results()[85].setName("Face_7_86")
Face_7.results()[86].setName("Face_7_87")
Face_7.results()[87].setName("Face_7_88")
Face_7.results()[88].setName("Face_7_89")
Face_7.results()[89].setName("Face_7_90")
Face_7.results()[90].setName("Face_7_91")
Face_7.results()[91].setName("Face_7_92")
Face_7.results()[92].setName("Face_7_93")
Face_7.results()[93].setName("Face_7_94")
Face_7.results()[94].setName("Face_7_95")
Face_7.results()[95].setName("Face_7_96")
Face_7.results()[96].setName("Face_7_97")
Face_7.results()[97].setName("Face_7_98")
Face_7.results()[98].setName("Face_7_99")
Face_7.results()[99].setName("Face_7_100")
Face_7.results()[100].setName("Face_7_101")
Face_7.results()[101].setName("Face_7_102")
Face_7.results()[102].setName("Face_7_103")
Face_7.results()[103].setName("Face_7_104")
Face_7.results()[104].setName("Face_7_105")
Face_7.results()[105].setName("Face_7_106")
Face_7.results()[106].setName("Face_7_107")
Face_7.results()[107].setName("Face_7_108")
Face_7.results()[108].setName("Face_7_109")
Face_7.results()[109].setName("Face_7_110")
Face_7.results()[110].setName("Face_7_111")
Face_7.results()[111].setName("Face_7_112")
Face_7.results()[112].setName("Face_7_113")
Face_7.results()[113].setName("Face_7_114")
Face_7.results()[114].setName("Face_7_115")
Face_7.results()[115].setName("Face_7_116")
Face_7.results()[116].setName("Face_7_117")
Face_7.results()[117].setName("Face_7_118")
Face_7.results()[118].setName("Face_7_119")
Face_7.results()[119].setName("Face_7_120")
Face_7.results()[120].setName("Face_7_121")
Face_7.results()[121].setName("Face_7_122")
Face_7.results()[122].setName("Face_7_123")
Face_7.results()[123].setName("Face_7_124")
Face_7.results()[124].setName("Face_7_125")
Face_7.results()[125].setName("Face_7_126")
Face_7.results()[126].setName("Face_7_127")
Face_7.results()[127].setName("Face_7_128")
Face_7.results()[128].setName("Face_7_129")
Face_7.results()[129].setName("Face_7_130")
Face_7.results()[130].setName("Face_7_131")
Face_7.results()[131].setName("Face_7_132")
Face_7.results()[132].setName("Face_7_133")
Face_7.results()[133].setName("Face_7_134")

### Create Sketch
Sketch_8 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_283 = Sketch_8.addLine(0.6734999999999807, -1.5128, 0.6735, -1.5033)

### Create SketchLine
SketchLine_284 = Sketch_8.addLine(0.674, -1.5028, 0.6779999999999999, -1.5028)

### Create SketchLine
SketchLine_285 = Sketch_8.addLine(0.6779999999999999, -1.5028, 0.6779999999999999, -1.4998)

### Create SketchLine
SketchLine_286 = Sketch_8.addLine(0.6779999999999999, -1.4998, 0.6799999999999999, -1.4998)

### Create SketchLine
SketchLine_287 = Sketch_8.addLine(0.6799999999999999, -1.4998, 0.6799999999999999, -1.499)

### Create SketchLine
SketchLine_288 = Sketch_8.addLine(0.6809999999999999, -1.498, 0.6909999999999241, -1.498)

### Create SketchLine
SketchLine_289 = Sketch_8.addLine(0.6929999999999241, -1.5, 0.6929999999999999, -1.5028)

### Create SketchLine
SketchLine_290 = Sketch_8.addLine(0.6929999999999999, -1.5028, 0.6970000000000001, -1.502800000000055)

### Create SketchLine
SketchLine_291 = Sketch_8.addLine(0.6975, -1.503300000000055, 0.6974999999999791, -1.5128)

### Create SketchLine
SketchLine_292 = Sketch_8.addLine(0.6874999999999791, -1.5228, 0.6834999999999807, -1.5228)
Sketch_8.setCoincident(SketchLine_284.endPoint(), SketchLine_285.startPoint())
Sketch_8.setCoincident(SketchLine_285.endPoint(), SketchLine_286.startPoint())
Sketch_8.setCoincident(SketchLine_286.endPoint(), SketchLine_287.startPoint())
Sketch_8.setCoincident(SketchLine_289.endPoint(), SketchLine_290.startPoint())

### Create SketchPoint
SketchPoint_19 = Sketch_8.addPoint(0.6735, -1.5128)
SketchPoint_19.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_19.result(), SketchLine_283.result())
Sketch_8.setFixed(SketchPoint_19.result())

### Create SketchPoint
SketchPoint_20 = Sketch_8.addPoint(0.67575, -1.5028)
SketchPoint_20.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_20.result(), SketchLine_284.result())
Sketch_8.setFixed(SketchPoint_20.result())

### Create SketchPoint
SketchPoint_21 = Sketch_8.addPoint(0.6779999999999999, -1.5013)
SketchPoint_21.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_21.result(), SketchLine_285.result())
Sketch_8.setFixed(SketchPoint_21.result())

### Create SketchPoint
SketchPoint_22 = Sketch_8.addPoint(0.6789999999999999, -1.4998)
SketchPoint_22.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_22.result(), SketchLine_286.result())
Sketch_8.setFixed(SketchPoint_22.result())

### Create SketchPoint
SketchPoint_23 = Sketch_8.addPoint(0.6799999999999999, -1.4989)
SketchPoint_23.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_23.result(), SketchLine_287.result())
Sketch_8.setFixed(SketchPoint_23.result())

### Create SketchPoint
SketchPoint_24 = Sketch_8.addPoint(0.6864999999999999, -1.498)
SketchPoint_24.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_24.result(), SketchLine_288.result())
Sketch_8.setFixed(SketchPoint_24.result())

### Create SketchPoint
SketchPoint_25 = Sketch_8.addPoint(0.6929999999999999, -1.5004)
SketchPoint_25.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_25.result(), SketchLine_289.result())
Sketch_8.setFixed(SketchPoint_25.result())

### Create SketchPoint
SketchPoint_26 = Sketch_8.addPoint(0.6952499999999999, -1.5028)
SketchPoint_26.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_26.result(), SketchLine_290.result())
Sketch_8.setFixed(SketchPoint_26.result())

### Create SketchPoint
SketchPoint_27 = Sketch_8.addPoint(0.6975, -1.5128)
SketchPoint_27.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_27.result(), SketchLine_291.result())
Sketch_8.setFixed(SketchPoint_27.result())

### Create SketchPoint
SketchPoint_28 = Sketch_8.addPoint(0.6855, -1.5228)
SketchPoint_28.setAuxiliary(True)
Sketch_8.setCoincident(SketchPoint_28.result(), SketchLine_292.result())
Sketch_8.setFixed(SketchPoint_28.result())

### Create SketchArc
SketchArc_270 = Sketch_8.addArc(0.674, -1.5033, 0.674, -1.5028, 0.6735, -1.5033, False)
Sketch_8.setCoincident(SketchArc_270.startPoint(), SketchLine_284.startPoint())
Sketch_8.setCoincident(SketchArc_270.endPoint(), SketchLine_283.endPoint())
Sketch_8.setTangent(SketchArc_270.results()[1], SketchLine_283.result())
Sketch_8.setTangent(SketchArc_270.results()[1], SketchLine_284.result())
Sketch_8.setRadius(SketchArc_270.results()[1], 0.0005)

### Create SketchArc
SketchArc_271 = Sketch_8.addArc(0.6970000000000001, -1.503300000000055, 0.6975, -1.503300000000055, 0.6970000000000001, -1.502800000000055, False)
Sketch_8.setCoincident(SketchArc_271.startPoint(), SketchLine_291.startPoint())
Sketch_8.setCoincident(SketchArc_271.endPoint(), SketchLine_290.endPoint())
Sketch_8.setTangent(SketchArc_271.results()[1], SketchLine_290.result())
Sketch_8.setTangent(SketchArc_271.results()[1], SketchLine_291.result())
Sketch_8.setRadius(SketchArc_271.results()[1], 0.0005)

### Create SketchArc
SketchArc_272 = Sketch_8.addArc(0.6834999999999807, -1.5128, 0.6734999999999807, -1.5128, 0.6834999999999807, -1.5228, False)
Sketch_8.setCoincident(SketchArc_272.startPoint(), SketchLine_283.startPoint())
Sketch_8.setCoincident(SketchArc_272.endPoint(), SketchLine_292.endPoint())
Sketch_8.setTangent(SketchArc_272.results()[1], SketchLine_292.result())
Sketch_8.setTangent(SketchArc_272.results()[1], SketchLine_283.result())
Sketch_8.setRadius(SketchArc_272.results()[1], 0.01)

### Create SketchArc
SketchArc_273 = Sketch_8.addArc(0.6874999999999791, -1.5128, 0.6874999999999791, -1.5228, 0.6974999999999791, -1.5128, False)
Sketch_8.setCoincident(SketchArc_273.startPoint(), SketchLine_292.startPoint())
Sketch_8.setCoincident(SketchArc_273.endPoint(), SketchLine_291.endPoint())
Sketch_8.setTangent(SketchArc_273.results()[1], SketchLine_292.result())
Sketch_8.setTangent(SketchArc_273.results()[1], SketchLine_291.result())
Sketch_8.setRadius(SketchArc_273.results()[1], 0.01)

### Create SketchArc
SketchArc_274 = Sketch_8.addArc(0.6909999999999241, -1.5, 0.6929999999999241, -1.5, 0.6909999999999241, -1.498, False)
Sketch_8.setCoincident(SketchArc_274.startPoint(), SketchLine_289.startPoint())
Sketch_8.setCoincident(SketchArc_274.endPoint(), SketchLine_288.endPoint())
Sketch_8.setTangent(SketchArc_274.results()[1], SketchLine_288.result())
Sketch_8.setTangent(SketchArc_274.results()[1], SketchLine_289.result())
Sketch_8.setRadius(SketchArc_274.results()[1], 0.002)

### Create SketchArc
SketchArc_275 = Sketch_8.addArc(0.6809999999999999, -1.499, 0.6809999999999999, -1.498, 0.6799999999999999, -1.499, False)
Sketch_8.setCoincident(SketchArc_275.startPoint(), SketchLine_288.startPoint())
Sketch_8.setCoincident(SketchArc_275.endPoint(), SketchLine_287.endPoint())
Sketch_8.setTangent(SketchArc_275.results()[1], SketchLine_288.result())
Sketch_8.setTangent(SketchArc_275.results()[1], SketchLine_287.result())
Sketch_8.setRadius(SketchArc_275.results()[1], 0.001)

### Create SketchCircle
SketchCircle_338 = Sketch_8.addCircle(0.07687499999999992, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_339 = Sketch_8.addCircle(0.556875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_340 = Sketch_8.addCircle(0.6243749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_341 = Sketch_8.addCircle(0.00187499999999996, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_342 = Sketch_8.addCircle(0.009375000000000022, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_343 = Sketch_8.addCircle(0.316875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_344 = Sketch_8.addCircle(0.631875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_345 = Sketch_8.addCircle(0.4818749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_346 = Sketch_8.addCircle(0.4893749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_347 = Sketch_8.addCircle(0.526875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_348 = Sketch_8.addCircle(0.384375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_349 = Sketch_8.addCircle(0.4218749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_350 = Sketch_8.addCircle(0.256875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_351 = Sketch_8.addCircle(0.271875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_352 = Sketch_8.addCircle(0.2793749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_353 = Sketch_8.addCircle(0.204375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_354 = Sketch_8.addCircle(0.2343749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_355 = Sketch_8.addCircle(0.249375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_356 = Sketch_8.addCircle(0.01687499999999997, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_357 = Sketch_8.addCircle(0.09187499999999993, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_358 = Sketch_8.addCircle(0.1068749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_359 = Sketch_8.addCircle(0.114375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_360 = Sketch_8.addCircle(0.1293749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_361 = Sketch_8.addCircle(0.444375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_362 = Sketch_8.addCircle(0.4668749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_363 = Sketch_8.addCircle(0.4743749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_364 = Sketch_8.addCircle(0.459375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_365 = Sketch_8.addCircle(0.3993749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_366 = Sketch_8.addCircle(0.661875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_367 = Sketch_8.addCircle(0.264375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_368 = Sketch_8.addCircle(0.219375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_369 = Sketch_8.addCircle(0.241875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_370 = Sketch_8.addCircle(0.2868749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_371 = Sketch_8.addCircle(0.504375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_372 = Sketch_8.addCircle(0.3468749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_373 = Sketch_8.addCircle(0.324375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_374 = Sketch_8.addCircle(0.3543749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_375 = Sketch_8.addCircle(0.331875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_376 = Sketch_8.addCircle(0.309375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_377 = Sketch_8.addCircle(0.436875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_378 = Sketch_8.addCircle(0.6543749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_379 = Sketch_8.addCircle(0.05437499999999995, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_380 = Sketch_8.addCircle(0.02437499999999992, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_381 = Sketch_8.addCircle(0.06187500000000001, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_382 = Sketch_8.addCircle(0.08437499999999998, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_383 = Sketch_8.addCircle(0.5343749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_384 = Sketch_8.addCircle(0.616875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_385 = Sketch_8.addCircle(0.5868749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_386 = Sketch_8.addCircle(0.6393749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_387 = Sketch_8.addCircle(0.03937499999999994, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_388 = Sketch_8.addCircle(0.136875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_389 = Sketch_8.addCircle(0.121875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_390 = Sketch_8.addCircle(0.3018749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_391 = Sketch_8.addCircle(0.1593749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_392 = Sketch_8.addCircle(0.2943749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_393 = Sketch_8.addCircle(0.646875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_394 = Sketch_8.addCircle(0.4518749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_395 = Sketch_8.addCircle(0.196875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_396 = Sketch_8.addCircle(0.1743749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_397 = Sketch_8.addCircle(0.166875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_398 = Sketch_8.addCircle(0.5418749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_399 = Sketch_8.addCircle(0.3618749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_400 = Sketch_8.addCircle(0.4143749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_401 = Sketch_8.addCircle(0.429375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_402 = Sketch_8.addCircle(0.376875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_403 = Sketch_8.addCircle(0.339375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_404 = Sketch_8.addCircle(0.391875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_405 = Sketch_8.addCircle(0.369375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_406 = Sketch_8.addCircle(0.4068749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_407 = Sketch_8.addCircle(0.181875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_408 = Sketch_8.addCircle(0.1443749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_409 = Sketch_8.addCircle(0.09937499999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_410 = Sketch_8.addCircle(0.189375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_411 = Sketch_8.addCircle(0.211875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_412 = Sketch_8.addCircle(0.2268749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_413 = Sketch_8.addCircle(0.151875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_414 = Sketch_8.addCircle(0.046875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_415 = Sketch_8.addCircle(0.06937499999999996, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_416 = Sketch_8.addCircle(0.03187499999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_417 = Sketch_8.addCircle(0.609375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_418 = Sketch_8.addCircle(0.5718749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_419 = Sketch_8.addCircle(0.594375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_420 = Sketch_8.addCircle(0.5193749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_421 = Sketch_8.addCircle(0.511875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_422 = Sketch_8.addCircle(0.496875, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_423 = Sketch_8.addCircle(0.564375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_424 = Sketch_8.addCircle(0.579375, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_425 = Sketch_8.addCircle(0.5493749999999999, -1.5028, 0.0001500000000000945)

### Create SketchCircle
SketchCircle_426 = Sketch_8.addCircle(0.6018749999999999, -1.5028, 0.0001500000000000945)
model.do()
Sketch_8.setName("Cathode")
Sketch_8.result().setName("Cathode")

### Create Face
Face_8 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "Cathode")])
Face_8.setName("Cathode")
Face_8.result().setName("Face_8_1")
Face_8.results()[1].setName("Face_8_2")
Face_8.results()[2].setName("Face_8_3")
Face_8.results()[3].setName("Face_8_4")
Face_8.results()[4].setName("Face_8_5")
Face_8.results()[5].setName("Face_8_6")
Face_8.results()[6].setName("Face_8_7")
Face_8.results()[7].setName("Face_8_8")
Face_8.results()[8].setName("Face_8_9")
Face_8.results()[9].setName("Face_8_10")
Face_8.results()[10].setName("Face_8_11")
Face_8.results()[11].setName("Face_8_12")
Face_8.results()[12].setName("Face_8_13")
Face_8.results()[13].setName("Face_8_14")
Face_8.results()[14].setName("Face_8_15")
Face_8.results()[15].setName("Face_8_16")
Face_8.results()[16].setName("Face_8_17")
Face_8.results()[17].setName("Face_8_18")
Face_8.results()[18].setName("Face_8_19")
Face_8.results()[19].setName("Face_8_20")
Face_8.results()[20].setName("Face_8_21")
Face_8.results()[21].setName("Face_8_22")
Face_8.results()[22].setName("Face_8_23")
Face_8.results()[23].setName("Face_8_24")
Face_8.results()[24].setName("Face_8_25")
Face_8.results()[25].setName("Face_8_26")
Face_8.results()[26].setName("Face_8_27")
Face_8.results()[27].setName("Face_8_28")
Face_8.results()[28].setName("Face_8_29")
Face_8.results()[29].setName("Face_8_30")
Face_8.results()[30].setName("Face_8_31")
Face_8.results()[31].setName("Face_8_32")
Face_8.results()[32].setName("Face_8_33")
Face_8.results()[33].setName("Face_8_34")
Face_8.results()[34].setName("Face_8_35")
Face_8.results()[35].setName("Face_8_36")
Face_8.results()[36].setName("Face_8_37")
Face_8.results()[37].setName("Face_8_38")
Face_8.results()[38].setName("Face_8_39")
Face_8.results()[39].setName("Face_8_40")
Face_8.results()[40].setName("Face_8_41")
Face_8.results()[41].setName("Face_8_42")
Face_8.results()[42].setName("Face_8_43")
Face_8.results()[43].setName("Face_8_44")
Face_8.results()[44].setName("Face_8_45")
Face_8.results()[45].setName("Face_8_46")
Face_8.results()[46].setName("Face_8_47")
Face_8.results()[47].setName("Face_8_48")
Face_8.results()[48].setName("Face_8_49")
Face_8.results()[49].setName("Face_8_50")
Face_8.results()[50].setName("Face_8_51")
Face_8.results()[51].setName("Face_8_52")
Face_8.results()[52].setName("Face_8_53")
Face_8.results()[53].setName("Face_8_54")
Face_8.results()[54].setName("Face_8_55")
Face_8.results()[55].setName("Face_8_56")
Face_8.results()[56].setName("Face_8_57")
Face_8.results()[57].setName("Face_8_58")
Face_8.results()[58].setName("Face_8_59")
Face_8.results()[59].setName("Face_8_60")
Face_8.results()[60].setName("Face_8_61")
Face_8.results()[61].setName("Face_8_62")
Face_8.results()[62].setName("Face_8_63")
Face_8.results()[63].setName("Face_8_64")
Face_8.results()[64].setName("Face_8_65")
Face_8.results()[65].setName("Face_8_66")
Face_8.results()[66].setName("Face_8_67")
Face_8.results()[67].setName("Face_8_68")
Face_8.results()[68].setName("Face_8_69")
Face_8.results()[69].setName("Face_8_70")
Face_8.results()[70].setName("Face_8_71")
Face_8.results()[71].setName("Face_8_72")
Face_8.results()[72].setName("Face_8_73")
Face_8.results()[73].setName("Face_8_74")
Face_8.results()[74].setName("Face_8_75")
Face_8.results()[75].setName("Face_8_76")
Face_8.results()[76].setName("Face_8_77")
Face_8.results()[77].setName("Face_8_78")
Face_8.results()[78].setName("Face_8_79")
Face_8.results()[79].setName("Face_8_80")
Face_8.results()[80].setName("Face_8_81")
Face_8.results()[81].setName("Face_8_82")
Face_8.results()[82].setName("Face_8_83")
Face_8.results()[83].setName("Face_8_84")
Face_8.results()[84].setName("Face_8_85")
Face_8.results()[85].setName("Face_8_86")
Face_8.results()[86].setName("Face_8_87")
Face_8.results()[87].setName("Face_8_88")
Face_8.results()[88].setName("Face_8_89")
Face_8.results()[89].setName("Face_8_90")

### Create Sketch
Sketch_9 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_293 = Sketch_9.addLine(0.667, 0.04159999999999988, 0.667, 0.0534)

### Create SketchLine
SketchLine_294 = Sketch_9.addLine(0.6686000000000001, 0.055, 0.6964000000000002, 0.055)

### Create SketchLine
SketchLine_295 = Sketch_9.addLine(0.6980000000000003, 0.0534, 0.6980000000000002, 0.04160000000000001)

### Create SketchLine
SketchLine_296 = Sketch_9.addLine(0.6964000000000001, 0.04, 0.6686000000000001, 0.03999999999999988)

### Create SketchPoint
SketchPoint_29 = Sketch_9.addPoint(0.667, 0.0475)
SketchPoint_29.setAuxiliary(True)
Sketch_9.setCoincident(SketchPoint_29.result(), SketchLine_293.result())
Sketch_9.setFixed(SketchPoint_29.result())

### Create SketchPoint
SketchPoint_30 = Sketch_9.addPoint(0.6825000000000001, 0.055)
SketchPoint_30.setAuxiliary(True)
Sketch_9.setCoincident(SketchPoint_30.result(), SketchLine_294.result())
Sketch_9.setFixed(SketchPoint_30.result())

### Create SketchPoint
SketchPoint_31 = Sketch_9.addPoint(0.6980000000000001, 0.0475)
SketchPoint_31.setAuxiliary(True)
Sketch_9.setCoincident(SketchPoint_31.result(), SketchLine_295.result())
Sketch_9.setFixed(SketchPoint_31.result())

### Create SketchPoint
SketchPoint_32 = Sketch_9.addPoint(0.6825000000000001, 0.04)
SketchPoint_32.setAuxiliary(True)
Sketch_9.setCoincident(SketchPoint_32.result(), SketchLine_296.result())
Sketch_9.setFixed(SketchPoint_32.result())

### Create SketchArc
SketchArc_276 = Sketch_9.addArc(0.6686000000000001, 0.04159999999999989, 0.667, 0.04159999999999988, 0.6686000000000001, 0.03999999999999988, False)
Sketch_9.setCoincident(SketchArc_276.startPoint(), SketchLine_293.startPoint())
Sketch_9.setCoincident(SketchArc_276.endPoint(), SketchLine_296.endPoint())
Sketch_9.setTangent(SketchArc_276.results()[1], SketchLine_293.result())
Sketch_9.setTangent(SketchArc_276.results()[1], SketchLine_296.result())
Sketch_9.setRadius(SketchArc_276.results()[1], 0.0016)

### Create SketchArc
SketchArc_277 = Sketch_9.addArc(0.6686000000000001, 0.0534, 0.6686000000000001, 0.055, 0.667, 0.0534, False)
Sketch_9.setCoincident(SketchArc_277.startPoint(), SketchLine_294.startPoint())
Sketch_9.setCoincident(SketchArc_277.endPoint(), SketchLine_293.endPoint())
Sketch_9.setTangent(SketchArc_277.results()[1], SketchLine_294.result())
Sketch_9.setTangent(SketchArc_277.results()[1], SketchLine_293.result())
Sketch_9.setRadius(SketchArc_277.results()[1], 0.0016)

### Create SketchArc
SketchArc_278 = Sketch_9.addArc(0.6964000000000002, 0.0534, 0.6980000000000003, 0.0534, 0.6964000000000002, 0.055, False)
Sketch_9.setCoincident(SketchArc_278.startPoint(), SketchLine_295.startPoint())
Sketch_9.setCoincident(SketchArc_278.endPoint(), SketchLine_294.endPoint())
Sketch_9.setTangent(SketchArc_278.results()[1], SketchLine_294.result())
Sketch_9.setTangent(SketchArc_278.results()[1], SketchLine_295.result())
Sketch_9.setRadius(SketchArc_278.results()[1], 0.0016)

### Create SketchArc
SketchArc_279 = Sketch_9.addArc(0.6964000000000001, 0.0416, 0.6964000000000001, 0.04, 0.6980000000000002, 0.04160000000000001, False)
Sketch_9.setCoincident(SketchArc_279.startPoint(), SketchLine_296.startPoint())
Sketch_9.setCoincident(SketchArc_279.endPoint(), SketchLine_295.endPoint())
Sketch_9.setTangent(SketchArc_279.results()[1], SketchLine_296.result())
Sketch_9.setTangent(SketchArc_279.results()[1], SketchLine_295.result())
Sketch_9.setRadius(SketchArc_279.results()[1], 0.0016)
model.do()
Sketch_9.setName("TopScreen")
Sketch_9.result().setName("TopScreen")

### Create Face
Face_9 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "TopScreen")])
Face_9.setName("TopScreen")
Face_9.result().setName("Face_9_1")

### Create Sketch
Sketch_10 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_297 = Sketch_10.addLine(0.6725, -1.5575, 0.6725, -1.5505)

### Create SketchLine
SketchLine_298 = Sketch_10.addLine(0.6799999999999999, -1.543, 0.6899999999999894, -1.54299999999998)

### Create SketchLine
SketchLine_299 = Sketch_10.addLine(0.6974999999999895, -1.550499999999984, 0.6974999999999903, -1.5575)

### Create SketchLine
SketchLine_300 = Sketch_10.addLine(0.6969999999999903, -1.558, 0.673, -1.558)

### Create SketchPoint
SketchPoint_33 = Sketch_10.addPoint(0.6725, -1.5505)
SketchPoint_33.setAuxiliary(True)
Sketch_10.setCoincident(SketchPoint_33.result(), SketchLine_297.result())
Sketch_10.setFixed(SketchPoint_33.result())

### Create SketchPoint
SketchPoint_34 = Sketch_10.addPoint(0.6850000000000001, -1.543)
SketchPoint_34.setAuxiliary(True)
Sketch_10.setCoincident(SketchPoint_34.result(), SketchLine_298.result())
Sketch_10.setFixed(SketchPoint_34.result())

### Create SketchPoint
SketchPoint_35 = Sketch_10.addPoint(0.6975, -1.5505)
SketchPoint_35.setAuxiliary(True)
Sketch_10.setCoincident(SketchPoint_35.result(), SketchLine_299.result())
Sketch_10.setFixed(SketchPoint_35.result())

### Create SketchPoint
SketchPoint_36 = Sketch_10.addPoint(0.6850000000000001, -1.558)
SketchPoint_36.setAuxiliary(True)
Sketch_10.setCoincident(SketchPoint_36.result(), SketchLine_300.result())
Sketch_10.setFixed(SketchPoint_36.result())

### Create SketchArc
SketchArc_280 = Sketch_10.addArc(0.6729999999999999, -1.5575, 0.6725, -1.5575, 0.673, -1.558, False)
Sketch_10.setCoincident(SketchArc_280.startPoint(), SketchLine_297.startPoint())
Sketch_10.setCoincident(SketchArc_280.endPoint(), SketchLine_300.endPoint())
Sketch_10.setTangent(SketchArc_280.results()[1], SketchLine_300.result())
Sketch_10.setTangent(SketchArc_280.results()[1], SketchLine_297.result())
Sketch_10.setRadius(SketchArc_280.results()[1], 0.0005)

### Create SketchArc
SketchArc_281 = Sketch_10.addArc(0.6969999999999903, -1.5575, 0.6969999999999903, -1.558, 0.6974999999999903, -1.5575, False)
Sketch_10.setCoincident(SketchArc_281.startPoint(), SketchLine_300.startPoint())
Sketch_10.setCoincident(SketchArc_281.endPoint(), SketchLine_299.endPoint())
Sketch_10.setTangent(SketchArc_281.results()[1], SketchLine_300.result())
Sketch_10.setTangent(SketchArc_281.results()[1], SketchLine_299.result())
Sketch_10.setRadius(SketchArc_281.results()[1], 0.0005)

### Create SketchArc
SketchArc_282 = Sketch_10.addArc(0.6799999999999999, -1.5505, 0.6799999999999999, -1.543, 0.6725, -1.5505, False)
Sketch_10.setCoincident(SketchArc_282.startPoint(), SketchLine_298.startPoint())
Sketch_10.setCoincident(SketchArc_282.endPoint(), SketchLine_297.endPoint())
Sketch_10.setTangent(SketchArc_282.results()[1], SketchLine_297.result())
Sketch_10.setTangent(SketchArc_282.results()[1], SketchLine_298.result())
Sketch_10.setRadius(SketchArc_282.results()[1], 0.0075)

### Create SketchArc
SketchArc_283 = Sketch_10.addArc(0.6899999999999894, -1.55049999999998, 0.6974999999999895, -1.550499999999984, 0.6899999999999894, -1.54299999999998, False)
Sketch_10.setCoincident(SketchArc_283.startPoint(), SketchLine_299.startPoint())
Sketch_10.setCoincident(SketchArc_283.endPoint(), SketchLine_298.endPoint())
Sketch_10.setTangent(SketchArc_283.results()[1], SketchLine_298.result())
Sketch_10.setTangent(SketchArc_283.results()[1], SketchLine_299.result())
Sketch_10.setRadius(SketchArc_283.results()[1], 0.0075)
model.do()
Sketch_10.setName("BottomScreen")
Sketch_10.result().setName("BottomScreen")

### Create Face
Face_10 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "BottomScreen")])
Face_10.setName("BottomScreen")
Face_10.result().setName("Face_10_1")

### Create Sketch
Sketch_11 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_301 = Sketch_11.addLine(0.679, -0.0335, 0.679, -0.0265)

### Create SketchLine
SketchLine_302 = Sketch_11.addLine(0.6805, -0.025, 0.6985000000000001, -0.025)

### Create SketchLine
SketchLine_303 = Sketch_11.addLine(0.7000000000000001, -0.02650000000000001, 0.7000000000000005, -0.0335)

### Create SketchLine
SketchLine_304 = Sketch_11.addLine(0.6985000000000006, -0.035, 0.6805, -0.035)

### Create SketchPoint
SketchPoint_37 = Sketch_11.addPoint(0.679, -0.03)
SketchPoint_37.setAuxiliary(True)
Sketch_11.setCoincident(SketchPoint_37.result(), SketchLine_301.result())
Sketch_11.setFixed(SketchPoint_37.result())

### Create SketchPoint
SketchPoint_38 = Sketch_11.addPoint(0.6895, -0.025)
SketchPoint_38.setAuxiliary(True)
Sketch_11.setCoincident(SketchPoint_38.result(), SketchLine_302.result())
Sketch_11.setFixed(SketchPoint_38.result())

### Create SketchPoint
SketchPoint_39 = Sketch_11.addPoint(0.7000000000000001, -0.03)
SketchPoint_39.setAuxiliary(True)
Sketch_11.setCoincident(SketchPoint_39.result(), SketchLine_303.result())
Sketch_11.setFixed(SketchPoint_39.result())

### Create SketchPoint
SketchPoint_40 = Sketch_11.addPoint(0.6895, -0.035)
SketchPoint_40.setAuxiliary(True)
Sketch_11.setCoincident(SketchPoint_40.result(), SketchLine_304.result())
Sketch_11.setFixed(SketchPoint_40.result())

### Create SketchArc
SketchArc_284 = Sketch_11.addArc(0.6805, -0.0335, 0.679, -0.0335, 0.6805, -0.035, False)
Sketch_11.setCoincident(SketchArc_284.startPoint(), SketchLine_301.startPoint())
Sketch_11.setCoincident(SketchArc_284.endPoint(), SketchLine_304.endPoint())
Sketch_11.setTangent(SketchArc_284.results()[1], SketchLine_304.result())
Sketch_11.setTangent(SketchArc_284.results()[1], SketchLine_301.result())
Sketch_11.setRadius(SketchArc_284.results()[1], 0.0015)

### Create SketchArc
SketchArc_285 = Sketch_11.addArc(0.6805, -0.0265, 0.6805, -0.025, 0.679, -0.0265, False)
Sketch_11.setCoincident(SketchArc_285.startPoint(), SketchLine_302.startPoint())
Sketch_11.setCoincident(SketchArc_285.endPoint(), SketchLine_301.endPoint())
Sketch_11.setTangent(SketchArc_285.results()[1], SketchLine_302.result())
Sketch_11.setTangent(SketchArc_285.results()[1], SketchLine_301.result())
Sketch_11.setRadius(SketchArc_285.results()[1], 0.0015)

### Create SketchArc
SketchArc_286 = Sketch_11.addArc(0.6985000000000001, -0.02650000000000001, 0.7000000000000001, -0.02650000000000001, 0.6985000000000001, -0.025, False)
Sketch_11.setCoincident(SketchArc_286.startPoint(), SketchLine_303.startPoint())
Sketch_11.setCoincident(SketchArc_286.endPoint(), SketchLine_302.endPoint())
Sketch_11.setTangent(SketchArc_286.results()[1], SketchLine_303.result())
Sketch_11.setTangent(SketchArc_286.results()[1], SketchLine_302.result())
Sketch_11.setRadius(SketchArc_286.results()[1], 0.0015)

### Create SketchArc
SketchArc_287 = Sketch_11.addArc(0.6985000000000006, -0.0335, 0.6985000000000006, -0.035, 0.7000000000000005, -0.0335, False)
Sketch_11.setCoincident(SketchArc_287.startPoint(), SketchLine_304.startPoint())
Sketch_11.setCoincident(SketchArc_287.endPoint(), SketchLine_303.endPoint())
Sketch_11.setTangent(SketchArc_287.results()[1], SketchLine_303.result())
Sketch_11.setTangent(SketchArc_287.results()[1], SketchLine_304.result())
Sketch_11.setRadius(SketchArc_287.results()[1], 0.0015)
model.do()
Sketch_11.setName("CopperRing")
Sketch_11.result().setName("CopperRing")

### Create Face
Face_11 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "CopperRing")])
Face_11.setName("CopperRing")
Face_11.result().setName("Face_11_1")

### Create Sketch
Sketch_12 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_305 = Sketch_12.addLine(0.61, -1.5629, 0.6860000000000001, -1.5629)

### Create SketchLine
SketchLine_306 = Sketch_12.addLine(0.6860000000000001, -1.5629, 0.6860000000000001, -1.5764)

### Create SketchArc
SketchArc_288 = Sketch_12.addArc(0.648, -1.5764, 0.6860000000000001, -1.5764, 0.67465, -1.603488327744621, True)

### Create SketchLine
SketchLine_307 = Sketch_12.addLine(0.67465, -1.603488327744621, 0.67465, -1.6769)

### Create SketchLine
SketchLine_308 = Sketch_12.addLine(0.67465, -1.6769, 0.6639, -1.6769)

### Create SketchLine
SketchLine_309 = Sketch_12.addLine(0.6639, -1.6769, 0.6639, -1.68645)

### Create SketchLine
SketchLine_310 = Sketch_12.addLine(0.6639, -1.68645, 0.6321, -1.68645)

### Create SketchLine
SketchLine_311 = Sketch_12.addLine(0.6321, -1.68645, 0.6321, -1.6769)

### Create SketchLine
SketchLine_312 = Sketch_12.addLine(0.6321, -1.6769, 0.6213500000000001, -1.6769)

### Create SketchLine
SketchLine_313 = Sketch_12.addLine(0.6213500000000001, -1.6769, 0.6213500000000001, -1.603488327744621)

### Create SketchArc
SketchArc_289 = Sketch_12.addArc(0.648, -1.5764, 0.6213500000000001, -1.603488327744621, 0.61, -1.5764, True)

### Create SketchLine
SketchLine_314 = Sketch_12.addLine(0.61, -1.5764, 0.61, -1.5629)
Sketch_12.setCoincident(SketchLine_305.endPoint(), SketchLine_306.startPoint())
Sketch_12.setCoincident(SketchLine_306.endPoint(), SketchArc_288.startPoint())
Sketch_12.setCoincident(SketchArc_288.endPoint(), SketchLine_307.startPoint())
Sketch_12.setCoincident(SketchLine_307.endPoint(), SketchLine_308.startPoint())
Sketch_12.setCoincident(SketchLine_308.endPoint(), SketchLine_309.startPoint())
Sketch_12.setCoincident(SketchLine_309.endPoint(), SketchLine_310.startPoint())
Sketch_12.setCoincident(SketchLine_310.endPoint(), SketchLine_311.startPoint())
Sketch_12.setCoincident(SketchLine_311.endPoint(), SketchLine_312.startPoint())
Sketch_12.setCoincident(SketchLine_312.endPoint(), SketchLine_313.startPoint())
Sketch_12.setCoincident(SketchLine_313.endPoint(), SketchArc_289.startPoint())
Sketch_12.setCoincident(SketchArc_289.endPoint(), SketchLine_314.startPoint())
Sketch_12.setCoincident(SketchLine_314.endPoint(), SketchLine_305.startPoint())

### Create SketchPoint
SketchPoint_41 = Sketch_12.addPoint(0.648, -1.5629)

### Create SketchPoint
SketchPoint_42 = Sketch_12.addPoint(0.5670000000000001, -1.5629)

### Create SketchMultiTranslation
SketchMultiTranslation_5_objects = [SketchLine_305.result(), SketchLine_306.result(), SketchArc_288.results()[1], SketchLine_307.result(), SketchLine_308.result(), SketchLine_309.result(), SketchLine_310.result(), SketchLine_311.result(), SketchLine_312.result(), SketchLine_313.result(), SketchArc_289.results()[1], SketchLine_314.result()]
SketchMultiTranslation_5 = Sketch_12.addTranslation(SketchMultiTranslation_5_objects, SketchPoint_41.coordinates(), SketchPoint_42.coordinates(), 8)
[SketchLine_315, SketchLine_316, SketchLine_317, SketchLine_318, SketchLine_319, SketchLine_320, SketchLine_321, SketchLine_322, SketchLine_323, SketchLine_324, SketchLine_325, SketchLine_326, SketchLine_327, SketchLine_328, SketchArc_290, SketchArc_291, SketchArc_292, SketchArc_293, SketchArc_294, SketchArc_295, SketchArc_296, SketchLine_329, SketchLine_330, SketchLine_331, SketchLine_332, SketchLine_333, SketchLine_334, SketchLine_335, SketchLine_336, SketchLine_337, SketchLine_338, SketchLine_339, SketchLine_340, SketchLine_341, SketchLine_342, SketchLine_343, SketchLine_344, SketchLine_345, SketchLine_346, SketchLine_347, SketchLine_348, SketchLine_349, SketchLine_350, SketchLine_351, SketchLine_352, SketchLine_353, SketchLine_354, SketchLine_355, SketchLine_356, SketchLine_357, SketchLine_358, SketchLine_359, SketchLine_360, SketchLine_361, SketchLine_362, SketchLine_363, SketchLine_364, SketchLine_365, SketchLine_366, SketchLine_367, SketchLine_368, SketchLine_369, SketchLine_370, SketchLine_371, SketchLine_372, SketchLine_373, SketchLine_374, SketchLine_375, SketchLine_376, SketchLine_377, SketchArc_297, SketchArc_298, SketchArc_299, SketchArc_300, SketchArc_301, SketchArc_302, SketchArc_303, SketchLine_378, SketchLine_379, SketchLine_380, SketchLine_381, SketchLine_382, SketchLine_383, SketchLine_384] = SketchMultiTranslation_5.translatedList()

### Create SketchPoint
SketchPoint_43 = Sketch_12.addPoint(0, -1.5629)

### Create SketchPoint
SketchPoint_44 = Sketch_12.addPoint(0, -1.68645)

### Create SketchLine
SketchLine_385 = Sketch_12.addLine(0, -1.68645, 0, -1.5629)

### Create SketchLine
SketchLine_386 = Sketch_12.addLine(0, -1.5629, 0.038, -1.5629)

### Create SketchLine
SketchLine_387 = Sketch_12.addLine(0.038, -1.5629, 0.038, -1.5764)

### Create SketchArc
SketchArc_304 = Sketch_12.addArc(0, -1.5764, 0.038, -1.5764, 0.02665, -1.603488327744621, True)

### Create SketchLine
SketchLine_388 = Sketch_12.addLine(0.02665, -1.603488327744621, 0.02665, -1.6769)

### Create SketchLine
SketchLine_389 = Sketch_12.addLine(0.02665, -1.6769, 0.0159, -1.6769)

### Create SketchLine
SketchLine_390 = Sketch_12.addLine(0.0159, -1.6769, 0.0159, -1.68645)

### Create SketchLine
SketchLine_391 = Sketch_12.addLine(0.0159, -1.68645, 0, -1.68645)
Sketch_12.setCoincident(SketchLine_385.endPoint(), SketchLine_386.startPoint())
Sketch_12.setCoincident(SketchLine_386.endPoint(), SketchLine_387.startPoint())
Sketch_12.setCoincident(SketchLine_387.endPoint(), SketchArc_304.startPoint())
Sketch_12.setCoincident(SketchArc_304.endPoint(), SketchLine_388.startPoint())
Sketch_12.setCoincident(SketchLine_388.endPoint(), SketchLine_389.startPoint())
Sketch_12.setCoincident(SketchLine_389.endPoint(), SketchLine_390.startPoint())
Sketch_12.setCoincident(SketchLine_390.endPoint(), SketchLine_391.startPoint())
Sketch_12.setCoincident(SketchLine_391.endPoint(), SketchLine_385.startPoint())
model.do()
Sketch_12.setName("BottomPMT")
Sketch_12.result().setName("BottomPMT")

### Create Face
Face_12 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "BottomPMT")])
Face_12.setName("BottomPMT")
Face_12.result().setName("Face_12_1")
Face_12.results()[1].setName("Face_12_2")
Face_12.results()[2].setName("Face_12_3")
Face_12.results()[3].setName("Face_12_4")
Face_12.results()[4].setName("Face_12_5")
Face_12.results()[5].setName("Face_12_6")
Face_12.results()[6].setName("Face_12_7")
Face_12.results()[7].setName("Face_12_8")
Face_12.results()[8].setName("Face_12_9")

### Create Sketch
Sketch_13 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_392 = Sketch_13.addLine(0.6321, 0.227438327744621, 0.6639, 0.227438327744621)

### Create SketchLine
SketchLine_393 = Sketch_13.addLine(0.6639, 0.227438327744621, 0.6639, 0.217888327744621)

### Create SketchLine
SketchLine_394 = Sketch_13.addLine(0.6639, 0.217888327744621, 0.67465, 0.217888327744621)

### Create SketchLine
SketchLine_395 = Sketch_13.addLine(0.67465, 0.217888327744621, 0.67465, 0.117388327744621)

### Create SketchArc
SketchArc_305 = Sketch_13.addArc(0.648, 0.09029999999999999, 0.67465, 0.117388327744621, 0.6860000000000001, 0.09029999999999999, True)

### Create SketchLine
SketchLine_396 = Sketch_13.addLine(0.6860000000000001, 0.09029999999999999, 0.6860000000000001, 0.07679999999999999)

### Create SketchLine
SketchLine_397 = Sketch_13.addLine(0.6860000000000001, 0.07679999999999999, 0.61, 0.07679999999999999)

### Create SketchLine
SketchLine_398 = Sketch_13.addLine(0.61, 0.07679999999999999, 0.61, 0.09029999999999999)

### Create SketchArc
SketchArc_306 = Sketch_13.addArc(0.648, 0.09029999999999999, 0.61, 0.09029999999999999, 0.6213500000000001, 0.117388327744621, True)

### Create SketchLine
SketchLine_399 = Sketch_13.addLine(0.6213500000000001, 0.117388327744621, 0.6213500000000001, 0.217888327744621)

### Create SketchLine
SketchLine_400 = Sketch_13.addLine(0.6213500000000001, 0.217888327744621, 0.6321, 0.217888327744621)

### Create SketchLine
SketchLine_401 = Sketch_13.addLine(0.6321, 0.217888327744621, 0.6321, 0.227438327744621)
Sketch_13.setCoincident(SketchLine_397.endPoint(), SketchLine_398.startPoint())
Sketch_13.setCoincident(SketchLine_398.endPoint(), SketchArc_306.startPoint())
Sketch_13.setCoincident(SketchArc_306.endPoint(), SketchLine_399.startPoint())
Sketch_13.setCoincident(SketchLine_399.endPoint(), SketchLine_400.startPoint())
Sketch_13.setCoincident(SketchLine_400.endPoint(), SketchLine_401.startPoint())
Sketch_13.setCoincident(SketchLine_401.endPoint(), SketchLine_392.startPoint())
Sketch_13.setCoincident(SketchLine_392.endPoint(), SketchLine_393.startPoint())
Sketch_13.setCoincident(SketchLine_393.endPoint(), SketchLine_394.startPoint())
Sketch_13.setCoincident(SketchLine_394.endPoint(), SketchLine_395.startPoint())
Sketch_13.setCoincident(SketchLine_395.endPoint(), SketchArc_305.startPoint())
Sketch_13.setCoincident(SketchArc_305.endPoint(), SketchLine_396.startPoint())
Sketch_13.setCoincident(SketchLine_396.endPoint(), SketchLine_397.startPoint())

### Create SketchPoint
SketchPoint_45 = Sketch_13.addPoint(0.648, 0.07679999999999999)

### Create SketchPoint
SketchPoint_46 = Sketch_13.addPoint(0.5670000000000001, 0.07679999999999999)

### Create SketchMultiTranslation
SketchMultiTranslation_6_objects = [SketchLine_392.result(), SketchLine_393.result(), SketchLine_394.result(), SketchLine_395.result(), SketchArc_305.results()[1], SketchLine_396.result(), SketchLine_397.result(), SketchLine_398.result(), SketchArc_306.results()[1], SketchLine_399.result(), SketchLine_400.result(), SketchLine_401.result()]
SketchMultiTranslation_6 = Sketch_13.addTranslation(SketchMultiTranslation_6_objects, SketchPoint_45.coordinates(), SketchPoint_46.coordinates(), 8)
[SketchLine_402, SketchLine_403, SketchLine_404, SketchLine_405, SketchLine_406, SketchLine_407, SketchLine_408, SketchLine_409, SketchLine_410, SketchLine_411, SketchLine_412, SketchLine_413, SketchLine_414, SketchLine_415, SketchLine_416, SketchLine_417, SketchLine_418, SketchLine_419, SketchLine_420, SketchLine_421, SketchLine_422, SketchLine_423, SketchLine_424, SketchLine_425, SketchLine_426, SketchLine_427, SketchLine_428, SketchLine_429, SketchArc_307, SketchArc_308, SketchArc_309, SketchArc_310, SketchArc_311, SketchArc_312, SketchArc_313, SketchLine_430, SketchLine_431, SketchLine_432, SketchLine_433, SketchLine_434, SketchLine_435, SketchLine_436, SketchLine_437, SketchLine_438, SketchLine_439, SketchLine_440, SketchLine_441, SketchLine_442, SketchLine_443, SketchLine_444, SketchLine_445, SketchLine_446, SketchLine_447, SketchLine_448, SketchLine_449, SketchLine_450, SketchArc_314, SketchArc_315, SketchArc_316, SketchArc_317, SketchArc_318, SketchArc_319, SketchArc_320, SketchLine_451, SketchLine_452, SketchLine_453, SketchLine_454, SketchLine_455, SketchLine_456, SketchLine_457, SketchLine_458, SketchLine_459, SketchLine_460, SketchLine_461, SketchLine_462, SketchLine_463, SketchLine_464, SketchLine_465, SketchLine_466, SketchLine_467, SketchLine_468, SketchLine_469, SketchLine_470, SketchLine_471] = SketchMultiTranslation_6.translatedList()

### Create SketchLine
SketchLine_472 = Sketch_13.addLine(0, 0.2274383277446209, 0, 0.07679999999999999)

### Create SketchLine
SketchLine_473 = Sketch_13.addLine(0, 0.07679999999999999, 0.038, 0.07679999999999999)

### Create SketchLine
SketchLine_474 = Sketch_13.addLine(0.038, 0.07679999999999999, 0.038, 0.09029999999999999)

### Create SketchArc
SketchArc_321 = Sketch_13.addArc(0, 0.09029999999999999, 0.02665, 0.1173883277446209, 0.038, 0.09029999999999999, True)

### Create SketchLine
SketchLine_475 = Sketch_13.addLine(0.02665, 0.1173883277446209, 0.02665, 0.2178883277446209)

### Create SketchLine
SketchLine_476 = Sketch_13.addLine(0.02665, 0.2178883277446209, 0.0159, 0.2178883277446209)

### Create SketchLine
SketchLine_477 = Sketch_13.addLine(0.0159, 0.2178883277446209, 0.0159, 0.2274383277446209)

### Create SketchLine
SketchLine_478 = Sketch_13.addLine(0.0159, 0.2274383277446209, 0, 0.2274383277446209)
Sketch_13.setCoincident(SketchLine_472.endPoint(), SketchLine_473.startPoint())
Sketch_13.setCoincident(SketchLine_473.endPoint(), SketchLine_474.startPoint())
Sketch_13.setCoincident(SketchLine_474.endPoint(), SketchArc_321.endPoint())
Sketch_13.setCoincident(SketchArc_321.startPoint(), SketchLine_475.startPoint())
Sketch_13.setCoincident(SketchLine_475.endPoint(), SketchLine_476.startPoint())
Sketch_13.setCoincident(SketchLine_476.endPoint(), SketchLine_477.startPoint())
Sketch_13.setCoincident(SketchLine_477.endPoint(), SketchLine_478.startPoint())
Sketch_13.setCoincident(SketchLine_478.endPoint(), SketchLine_472.startPoint())
model.do()
Sketch_13.setName("TopPMT")
Sketch_13.result().setName("TopPMT")

### Create Face
Face_13 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "TopPMT")])
Face_13.setName("TopPMT")
Face_13.result().setName("Face_13_1")
Face_13.results()[1].setName("Face_13_2")
Face_13.results()[2].setName("Face_13_3")
Face_13.results()[3].setName("Face_13_4")
Face_13.results()[4].setName("Face_13_5")
Face_13.results()[5].setName("Face_13_6")
Face_13.results()[6].setName("Face_13_7")
Face_13.results()[7].setName("Face_13_8")
Face_13.results()[8].setName("Face_13_9")

### Create Sketch
Sketch_14 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_479 = Sketch_14.addLine(0.6643, 0.0005, 0.6643, 0.0075)

### Create SketchLine
SketchLine_480 = Sketch_14.addLine(0.6643, 0.0075, 0.6763, 0.0075)

### Create SketchLine
SketchLine_481 = Sketch_14.addLine(0.6763, 0.0075, 0.6763, 0.008)

### Create SketchLine
SketchLine_482 = Sketch_14.addLine(0.6763, 0.008, 0.6965, 0.008)

### Create SketchLine
SketchLine_483 = Sketch_14.addLine(0.6965, 0.008, 0.6965, 0.006)

### Create SketchLine
SketchLine_484 = Sketch_14.addLine(0.6965, 0.006, 0.7005, 0.006)

### Create SketchLine
SketchLine_485 = Sketch_14.addLine(0.7005, 0.006, 0.7005, -0.02)

### Create SketchLine
SketchLine_486 = Sketch_14.addLine(0.7005, -0.02, 0.6985, -0.02)

### Create SketchLine
SketchLine_487 = Sketch_14.addLine(0.6985, -0.02, 0.6985, 3.469446951953614e-18)

### Create SketchLine
SketchLine_488 = Sketch_14.addLine(0.6985, 3.469446951953614e-18, 0.6763, 0)

### Create SketchLine
SketchLine_489 = Sketch_14.addLine(0.6763, 0, 0.6763, 0.0005)

### Create SketchLine
SketchLine_490 = Sketch_14.addLine(0.6763, 0.0005, 0.6643, 0.0005)
Sketch_14.setCoincident(SketchLine_479.endPoint(), SketchLine_480.startPoint())
Sketch_14.setCoincident(SketchLine_480.endPoint(), SketchLine_481.startPoint())
Sketch_14.setCoincident(SketchLine_481.endPoint(), SketchLine_482.startPoint())
Sketch_14.setCoincident(SketchLine_482.endPoint(), SketchLine_483.startPoint())
Sketch_14.setCoincident(SketchLine_483.endPoint(), SketchLine_484.startPoint())
Sketch_14.setCoincident(SketchLine_484.endPoint(), SketchLine_485.startPoint())
Sketch_14.setCoincident(SketchLine_485.endPoint(), SketchLine_486.startPoint())
Sketch_14.setCoincident(SketchLine_486.endPoint(), SketchLine_487.startPoint())
Sketch_14.setCoincident(SketchLine_487.endPoint(), SketchLine_488.startPoint())
Sketch_14.setCoincident(SketchLine_488.endPoint(), SketchLine_489.startPoint())
Sketch_14.setCoincident(SketchLine_489.endPoint(), SketchLine_490.startPoint())
Sketch_14.setCoincident(SketchLine_490.endPoint(), SketchLine_479.startPoint())
model.do()
Sketch_14.setName("GateInsulatingFrame")
Sketch_14.result().setName("GateInsulatingFrame")

### Create Face
Face_14 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "GateInsulatingFrame")])
Face_14.setName("GateInsulatingFrame")
Face_14.result().setName("Face_14_1")

### Create Sketch
Sketch_15 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_491 = Sketch_15.addLine(0.6643, 0.008699999999999999, 0.6643, 0.0395)

### Create SketchLine
SketchLine_492 = Sketch_15.addLine(0.6643, 0.0395, 0.6643, 0.0543)

### Create SketchLine
SketchLine_493 = Sketch_15.addLine(0.6643, 0.0543, 0.667, 0.0543)

### Create SketchLine
SketchLine_494 = Sketch_15.addLine(0.667, 0.0543, 0.667, 0.0395)

### Create SketchLine
SketchLine_495 = Sketch_15.addLine(0.667, 0.0395, 0.6763, 0.0395)

### Create SketchLine
SketchLine_496 = Sketch_15.addLine(0.6763, 0.0395, 0.6763, 0.04)

### Create SketchLine
SketchLine_497 = Sketch_15.addLine(0.6763, 0.04, 0.6965, 0.04)

### Create SketchLine
SketchLine_498 = Sketch_15.addLine(0.6965, 0.04, 0.6965, 0.038)

### Create SketchLine
SketchLine_499 = Sketch_15.addLine(0.6965, 0.038, 0.7005, 0.038)

### Create SketchLine
SketchLine_500 = Sketch_15.addLine(0.7005, 0.038, 0.7005, 0.006499999999999999)

### Create SketchLine
SketchLine_501 = Sketch_15.addLine(0.7005, 0.006499999999999999, 0.6985, 0.006499999999999999)

### Create SketchLine
SketchLine_502 = Sketch_15.addLine(0.6985, 0.006499999999999999, 0.6985, 0.032)

### Create SketchLine
SketchLine_503 = Sketch_15.addLine(0.6985, 0.032, 0.6763, 0.032)

### Create SketchLine
SketchLine_504 = Sketch_15.addLine(0.6763, 0.032, 0.6763, 0.0325)

### Create SketchLine
SketchLine_505 = Sketch_15.addLine(0.6763, 0.0325, 0.667, 0.0325)

### Create SketchLine
SketchLine_506 = Sketch_15.addLine(0.667, 0.0325, 0.667, 0.008699999999999999)

### Create SketchLine
SketchLine_507 = Sketch_15.addLine(0.667, 0.008699999999999999, 0.6643, 0.008699999999999999)
Sketch_15.setCoincident(SketchLine_491.endPoint(), SketchLine_492.startPoint())
Sketch_15.setCoincident(SketchLine_492.endPoint(), SketchLine_493.startPoint())
Sketch_15.setCoincident(SketchLine_493.endPoint(), SketchLine_494.startPoint())
Sketch_15.setCoincident(SketchLine_494.endPoint(), SketchLine_495.startPoint())
Sketch_15.setCoincident(SketchLine_495.endPoint(), SketchLine_496.startPoint())
Sketch_15.setCoincident(SketchLine_496.endPoint(), SketchLine_497.startPoint())
Sketch_15.setCoincident(SketchLine_497.endPoint(), SketchLine_498.startPoint())
Sketch_15.setCoincident(SketchLine_498.endPoint(), SketchLine_499.startPoint())
Sketch_15.setCoincident(SketchLine_499.endPoint(), SketchLine_500.startPoint())
Sketch_15.setCoincident(SketchLine_500.endPoint(), SketchLine_501.startPoint())
Sketch_15.setCoincident(SketchLine_501.endPoint(), SketchLine_502.startPoint())
Sketch_15.setCoincident(SketchLine_502.endPoint(), SketchLine_503.startPoint())
Sketch_15.setCoincident(SketchLine_503.endPoint(), SketchLine_504.startPoint())
Sketch_15.setCoincident(SketchLine_504.endPoint(), SketchLine_505.startPoint())
Sketch_15.setCoincident(SketchLine_505.endPoint(), SketchLine_506.startPoint())
Sketch_15.setCoincident(SketchLine_506.endPoint(), SketchLine_507.startPoint())
Sketch_15.setCoincident(SketchLine_507.endPoint(), SketchLine_491.startPoint())
model.do()
Sketch_15.setName("AnodeInsulatingFrame")
Sketch_15.result().setName("AnodeInsulatingFrame")

### Create Face
Face_15 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "AnodeInsulatingFrame")])
Face_15.setName("AnodeInsulatingFrame")
Face_15.result().setName("Face_15_1")

### Create Sketch
Sketch_16 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_508 = Sketch_16.addLine(0.6691400000000001, -1.4763, 0.6721400000000001, -1.4763)

### Create SketchLine
SketchLine_509 = Sketch_16.addLine(0.6721400000000001, -1.4763, 0.6721400000000001, -1.4953)

### Create SketchLine
SketchLine_510 = Sketch_16.addLine(0.6721400000000001, -1.4953, 0.69964, -1.4953)

### Create SketchLine
SketchLine_511 = Sketch_16.addLine(0.69964, -1.4953, 0.69964, -1.5083)

### Create SketchLine
SketchLine_512 = Sketch_16.addLine(0.69964, -1.5083, 0.69764, -1.5083)

### Create SketchLine
SketchLine_513 = Sketch_16.addLine(0.69764, -1.5083, 0.69764, -1.5028)

### Create SketchLine
SketchLine_514 = Sketch_16.addLine(0.69764, -1.5028, 0.69564, -1.5028)

### Create SketchLine
SketchLine_515 = Sketch_16.addLine(0.69564, -1.5028, 0.69564, -1.49794)

### Create SketchLine
SketchLine_516 = Sketch_16.addLine(0.69564, -1.49794, 0.67764, -1.49794)

### Create SketchLine
SketchLine_517 = Sketch_16.addLine(0.67764, -1.49794, 0.67764, -1.5028)

### Create SketchLine
SketchLine_518 = Sketch_16.addLine(0.67764, -1.5028, 0.6691400000000001, -1.5028)

### Create SketchLine
SketchLine_519 = Sketch_16.addLine(0.6691400000000001, -1.5028, 0.6691400000000001, -1.4763)
Sketch_16.setCoincident(SketchLine_508.endPoint(), SketchLine_509.startPoint())
Sketch_16.setCoincident(SketchLine_509.endPoint(), SketchLine_510.startPoint())
Sketch_16.setCoincident(SketchLine_510.endPoint(), SketchLine_511.startPoint())
Sketch_16.setCoincident(SketchLine_511.endPoint(), SketchLine_512.startPoint())
Sketch_16.setCoincident(SketchLine_512.endPoint(), SketchLine_513.startPoint())
Sketch_16.setCoincident(SketchLine_513.endPoint(), SketchLine_514.startPoint())
Sketch_16.setCoincident(SketchLine_514.endPoint(), SketchLine_515.startPoint())
Sketch_16.setCoincident(SketchLine_515.endPoint(), SketchLine_516.startPoint())
Sketch_16.setCoincident(SketchLine_516.endPoint(), SketchLine_517.startPoint())
Sketch_16.setCoincident(SketchLine_517.endPoint(), SketchLine_518.startPoint())
Sketch_16.setCoincident(SketchLine_518.endPoint(), SketchLine_519.startPoint())
Sketch_16.setCoincident(SketchLine_519.endPoint(), SketchLine_508.startPoint())
model.do()
Sketch_16.setName("CathodeInsulatingFrame")
Sketch_16.result().setName("CathodeInsulatingFrame")

### Create Face
Face_16 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "CathodeInsulatingFrame")])
Face_16.setName("CathodeInsulatingFrame")
Face_16.result().setName("Face_16_1")

### Create Sketch
Sketch_17 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_520 = Sketch_17.addLine(0.6643, 0.07389999999999999, 0.667, 0.07389999999999999)

### Create SketchLine
SketchLine_521 = Sketch_17.addLine(0.667, 0.07389999999999999, 0.667, 0.059)

### Create SketchLine
SketchLine_522 = Sketch_17.addLine(0.667, 0.059, 0.7005, 0.059)

### Create SketchLine
SketchLine_523 = Sketch_17.addLine(0.7005, 0.059, 0.7005, 0.03849999999999999)

### Create SketchLine
SketchLine_524 = Sketch_17.addLine(0.7005, 0.03849999999999999, 0.6985, 0.03849999999999999)

### Create SketchLine
SketchLine_525 = Sketch_17.addLine(0.6985, 0.03849999999999999, 0.6985, 0.05499999999999999)

### Create SketchLine
SketchLine_526 = Sketch_17.addLine(0.6985, 0.05499999999999999, 0.6763, 0.05499999999999999)

### Create SketchLine
SketchLine_527 = Sketch_17.addLine(0.6763, 0.05499999999999999, 0.6763, 0.05549999999999999)

### Create SketchLine
SketchLine_528 = Sketch_17.addLine(0.6763, 0.05549999999999999, 0.667, 0.05549999999999999)

### Create SketchLine
SketchLine_529 = Sketch_17.addLine(0.667, 0.05549999999999999, 0.6643, 0.05549999999999999)

### Create SketchLine
SketchLine_530 = Sketch_17.addLine(0.6643, 0.05549999999999999, 0.6643, 0.07389999999999999)
Sketch_17.setCoincident(SketchLine_520.endPoint(), SketchLine_521.startPoint())
Sketch_17.setCoincident(SketchLine_521.endPoint(), SketchLine_522.startPoint())
Sketch_17.setCoincident(SketchLine_522.endPoint(), SketchLine_523.startPoint())
Sketch_17.setCoincident(SketchLine_523.endPoint(), SketchLine_524.startPoint())
Sketch_17.setCoincident(SketchLine_524.endPoint(), SketchLine_525.startPoint())
Sketch_17.setCoincident(SketchLine_525.endPoint(), SketchLine_526.startPoint())
Sketch_17.setCoincident(SketchLine_526.endPoint(), SketchLine_527.startPoint())
Sketch_17.setCoincident(SketchLine_527.endPoint(), SketchLine_528.startPoint())
Sketch_17.setCoincident(SketchLine_528.endPoint(), SketchLine_529.startPoint())
Sketch_17.setCoincident(SketchLine_529.endPoint(), SketchLine_530.startPoint())
Sketch_17.setCoincident(SketchLine_530.endPoint(), SketchLine_520.startPoint())
model.do()
Sketch_17.setName("TopScreenInsulatingFrame")
Sketch_17.result().setName("TopScreenInsulatingFrame")

### Create Face
Face_17 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "TopScreenInsulatingFrame")])
Face_17.setName("TopScreenInsulatingFrame")
Face_17.result().setName("Face_17_1")

### Create Sketch
Sketch_18 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchPoint
SketchPoint_47 = Sketch_18.addPoint(0.664, -1.50528)

### Create SketchPoint
SketchPoint_48 = Sketch_18.addPoint(0.664, -1.55528)

### Create SketchArc
SketchArc_322 = Sketch_18.addArc(0.664, -1.50528, 0.664, -1.50478, 0.664, -1.50578, True)

### Create SketchPoint
SketchPoint_49 = Sketch_18.addPoint(0.6682, -1.50528)

### Create SketchPoint
SketchPoint_50 = Sketch_18.addPoint(0.6682, -1.55528)

### Create SketchArc
SketchArc_323 = Sketch_18.addArc(0.6682, -1.55528, 0.6682, -1.55578, 0.6682, -1.55478, True)

### Create SketchMultiTranslation
SketchMultiTranslation_7 = Sketch_18.addTranslation([SketchArc_322.results()[1]], SketchPoint_47.coordinates(), SketchPoint_48.coordinates(), 26, True)
[SketchArc_324, SketchArc_325, SketchArc_326, SketchArc_327, SketchArc_328, SketchArc_329, SketchArc_330, SketchArc_331, SketchArc_332, SketchArc_333, SketchArc_334, SketchArc_335, SketchArc_336, SketchArc_337, SketchArc_338, SketchArc_339, SketchArc_340, SketchArc_341, SketchArc_342, SketchArc_343, SketchArc_344, SketchArc_345, SketchArc_346, SketchArc_347, SketchArc_348] = SketchMultiTranslation_7.translatedList()

### Create SketchMultiTranslation
SketchMultiTranslation_8 = Sketch_18.addTranslation([SketchArc_323.results()[1]], SketchPoint_50.coordinates(), SketchPoint_49.coordinates(), 26, True)
[SketchArc_349, SketchArc_350, SketchArc_351, SketchArc_352, SketchArc_353, SketchArc_354, SketchArc_355, SketchArc_356, SketchArc_357, SketchArc_358, SketchArc_359, SketchArc_360, SketchArc_361, SketchArc_362, SketchArc_363, SketchArc_364, SketchArc_365, SketchArc_366, SketchArc_367, SketchArc_368, SketchArc_369, SketchArc_370, SketchArc_371, SketchArc_372, SketchArc_373] = SketchMultiTranslation_8.translatedList()

### Create SketchLine
SketchLine_531 = Sketch_18.addLine(0.664, -1.50578, 0.664, -1.50678)
Sketch_18.setCoincident(SketchArc_322.endPoint(), SketchLine_531.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_324).startPoint(), SketchLine_531.endPoint())

### Create SketchLine
SketchLine_532 = Sketch_18.addLine(0.664, -1.50778, 0.664, -1.50878)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_324).endPoint(), SketchLine_532.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_325).startPoint(), SketchLine_532.endPoint())

### Create SketchLine
SketchLine_533 = Sketch_18.addLine(0.664, -1.50978, 0.664, -1.51078)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_325).endPoint(), SketchLine_533.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_326).startPoint(), SketchLine_533.endPoint())

### Create SketchLine
SketchLine_534 = Sketch_18.addLine(0.664, -1.51178, 0.664, -1.51278)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_326).endPoint(), SketchLine_534.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_327).startPoint(), SketchLine_534.endPoint())

### Create SketchLine
SketchLine_535 = Sketch_18.addLine(0.664, -1.51378, 0.664, -1.51478)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_327).endPoint(), SketchLine_535.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_328).startPoint(), SketchLine_535.endPoint())

### Create SketchLine
SketchLine_536 = Sketch_18.addLine(0.664, -1.51578, 0.664, -1.51678)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_328).endPoint(), SketchLine_536.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_329).startPoint(), SketchLine_536.endPoint())

### Create SketchLine
SketchLine_537 = Sketch_18.addLine(0.664, -1.51778, 0.664, -1.51878)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_329).endPoint(), SketchLine_537.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_330).startPoint(), SketchLine_537.endPoint())

### Create SketchLine
SketchLine_538 = Sketch_18.addLine(0.664, -1.51978, 0.664, -1.52078)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_330).endPoint(), SketchLine_538.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_331).startPoint(), SketchLine_538.endPoint())

### Create SketchLine
SketchLine_539 = Sketch_18.addLine(0.664, -1.52178, 0.664, -1.52278)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_331).endPoint(), SketchLine_539.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_332).startPoint(), SketchLine_539.endPoint())

### Create SketchLine
SketchLine_540 = Sketch_18.addLine(0.664, -1.52378, 0.664, -1.52478)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_332).endPoint(), SketchLine_540.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_333).startPoint(), SketchLine_540.endPoint())

### Create SketchLine
SketchLine_541 = Sketch_18.addLine(0.664, -1.52578, 0.664, -1.52678)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_333).endPoint(), SketchLine_541.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_334).startPoint(), SketchLine_541.endPoint())

### Create SketchLine
SketchLine_542 = Sketch_18.addLine(0.664, -1.52778, 0.664, -1.52878)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_334).endPoint(), SketchLine_542.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_335).startPoint(), SketchLine_542.endPoint())

### Create SketchLine
SketchLine_543 = Sketch_18.addLine(0.664, -1.52978, 0.664, -1.53078)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_335).endPoint(), SketchLine_543.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_336).startPoint(), SketchLine_543.endPoint())

### Create SketchLine
SketchLine_544 = Sketch_18.addLine(0.664, -1.53178, 0.664, -1.53278)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_336).endPoint(), SketchLine_544.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_337).startPoint(), SketchLine_544.endPoint())

### Create SketchLine
SketchLine_545 = Sketch_18.addLine(0.664, -1.53378, 0.664, -1.53478)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_337).endPoint(), SketchLine_545.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_338).startPoint(), SketchLine_545.endPoint())

### Create SketchLine
SketchLine_546 = Sketch_18.addLine(0.664, -1.53578, 0.664, -1.53678)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_338).endPoint(), SketchLine_546.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_339).startPoint(), SketchLine_546.endPoint())

### Create SketchLine
SketchLine_547 = Sketch_18.addLine(0.664, -1.53778, 0.664, -1.53878)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_339).endPoint(), SketchLine_547.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_340).startPoint(), SketchLine_547.endPoint())

### Create SketchLine
SketchLine_548 = Sketch_18.addLine(0.664, -1.53978, 0.664, -1.54078)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_340).endPoint(), SketchLine_548.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_341).startPoint(), SketchLine_548.endPoint())

### Create SketchLine
SketchLine_549 = Sketch_18.addLine(0.664, -1.54178, 0.664, -1.54278)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_341).endPoint(), SketchLine_549.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_342).startPoint(), SketchLine_549.endPoint())

### Create SketchLine
SketchLine_550 = Sketch_18.addLine(0.664, -1.54378, 0.664, -1.54478)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_342).endPoint(), SketchLine_550.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_343).startPoint(), SketchLine_550.endPoint())

### Create SketchLine
SketchLine_551 = Sketch_18.addLine(0.664, -1.54578, 0.664, -1.54678)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_343).endPoint(), SketchLine_551.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_344).startPoint(), SketchLine_551.endPoint())

### Create SketchLine
SketchLine_552 = Sketch_18.addLine(0.664, -1.54778, 0.664, -1.54878)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_344).endPoint(), SketchLine_552.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_345).startPoint(), SketchLine_552.endPoint())

### Create SketchLine
SketchLine_553 = Sketch_18.addLine(0.664, -1.54978, 0.664, -1.55078)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_345).endPoint(), SketchLine_553.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_346).startPoint(), SketchLine_553.endPoint())

### Create SketchLine
SketchLine_554 = Sketch_18.addLine(0.664, -1.55178, 0.664, -1.55278)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_346).endPoint(), SketchLine_554.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_347).startPoint(), SketchLine_554.endPoint())

### Create SketchLine
SketchLine_555 = Sketch_18.addLine(0.664, -1.55378, 0.664, -1.55478)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_347).endPoint(), SketchLine_555.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_348).startPoint(), SketchLine_555.endPoint())

### Create SketchLine
SketchLine_556 = Sketch_18.addLine(0.6682, -1.55478, 0.6682, -1.55378)
Sketch_18.setCoincident(SketchArc_323.endPoint(), SketchLine_556.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_349).startPoint(), SketchLine_556.endPoint())

### Create SketchLine
SketchLine_557 = Sketch_18.addLine(0.6682, -1.55278, 0.6682, -1.55178)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_349).endPoint(), SketchLine_557.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_350).startPoint(), SketchLine_557.endPoint())

### Create SketchLine
SketchLine_558 = Sketch_18.addLine(0.6682, -1.55078, 0.6682, -1.54978)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_350).endPoint(), SketchLine_558.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_351).startPoint(), SketchLine_558.endPoint())

### Create SketchLine
SketchLine_559 = Sketch_18.addLine(0.6682, -1.54878, 0.6682, -1.54778)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_351).endPoint(), SketchLine_559.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_352).startPoint(), SketchLine_559.endPoint())

### Create SketchLine
SketchLine_560 = Sketch_18.addLine(0.6682, -1.54678, 0.6682, -1.54578)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_352).endPoint(), SketchLine_560.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_353).startPoint(), SketchLine_560.endPoint())

### Create SketchLine
SketchLine_561 = Sketch_18.addLine(0.6682, -1.54478, 0.6682, -1.54378)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_353).endPoint(), SketchLine_561.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_354).startPoint(), SketchLine_561.endPoint())

### Create SketchLine
SketchLine_562 = Sketch_18.addLine(0.6682, -1.54278, 0.6682, -1.54178)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_354).endPoint(), SketchLine_562.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_355).startPoint(), SketchLine_562.endPoint())

### Create SketchLine
SketchLine_563 = Sketch_18.addLine(0.6682, -1.54078, 0.6682, -1.53978)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_355).endPoint(), SketchLine_563.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_356).startPoint(), SketchLine_563.endPoint())

### Create SketchLine
SketchLine_564 = Sketch_18.addLine(0.6682, -1.53878, 0.6682, -1.53778)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_356).endPoint(), SketchLine_564.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_357).startPoint(), SketchLine_564.endPoint())

### Create SketchLine
SketchLine_565 = Sketch_18.addLine(0.6682, -1.53678, 0.6682, -1.53578)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_357).endPoint(), SketchLine_565.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_358).startPoint(), SketchLine_565.endPoint())

### Create SketchLine
SketchLine_566 = Sketch_18.addLine(0.6682, -1.53478, 0.6682, -1.53378)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_358).endPoint(), SketchLine_566.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_359).startPoint(), SketchLine_566.endPoint())

### Create SketchLine
SketchLine_567 = Sketch_18.addLine(0.6682, -1.53278, 0.6682, -1.53178)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_359).endPoint(), SketchLine_567.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_360).startPoint(), SketchLine_567.endPoint())

### Create SketchLine
SketchLine_568 = Sketch_18.addLine(0.6682, -1.53078, 0.6682, -1.52978)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_360).endPoint(), SketchLine_568.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_361).startPoint(), SketchLine_568.endPoint())

### Create SketchLine
SketchLine_569 = Sketch_18.addLine(0.6682, -1.52878, 0.6682, -1.52778)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_361).endPoint(), SketchLine_569.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_362).startPoint(), SketchLine_569.endPoint())

### Create SketchLine
SketchLine_570 = Sketch_18.addLine(0.6682, -1.52678, 0.6682, -1.52578)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_362).endPoint(), SketchLine_570.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_363).startPoint(), SketchLine_570.endPoint())

### Create SketchLine
SketchLine_571 = Sketch_18.addLine(0.6682, -1.52478, 0.6682, -1.52378)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_363).endPoint(), SketchLine_571.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_364).startPoint(), SketchLine_571.endPoint())

### Create SketchLine
SketchLine_572 = Sketch_18.addLine(0.6682, -1.52278, 0.6682, -1.52178)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_364).endPoint(), SketchLine_572.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_365).startPoint(), SketchLine_572.endPoint())

### Create SketchLine
SketchLine_573 = Sketch_18.addLine(0.6682, -1.52078, 0.6682, -1.51978)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_365).endPoint(), SketchLine_573.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_366).startPoint(), SketchLine_573.endPoint())

### Create SketchLine
SketchLine_574 = Sketch_18.addLine(0.6682, -1.51878, 0.6682, -1.51778)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_366).endPoint(), SketchLine_574.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_367).startPoint(), SketchLine_574.endPoint())

### Create SketchLine
SketchLine_575 = Sketch_18.addLine(0.6682, -1.51678, 0.6682, -1.51578)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_367).endPoint(), SketchLine_575.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_368).startPoint(), SketchLine_575.endPoint())

### Create SketchLine
SketchLine_576 = Sketch_18.addLine(0.6682, -1.51478, 0.6682, -1.51378)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_368).endPoint(), SketchLine_576.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_369).startPoint(), SketchLine_576.endPoint())

### Create SketchLine
SketchLine_577 = Sketch_18.addLine(0.6682, -1.51278, 0.6682, -1.51178)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_369).endPoint(), SketchLine_577.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_370).startPoint(), SketchLine_577.endPoint())

### Create SketchLine
SketchLine_578 = Sketch_18.addLine(0.6682, -1.51078, 0.6682, -1.50978)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_370).endPoint(), SketchLine_578.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_371).startPoint(), SketchLine_578.endPoint())

### Create SketchLine
SketchLine_579 = Sketch_18.addLine(0.6682, -1.50878, 0.6682, -1.50778)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_371).endPoint(), SketchLine_579.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_372).startPoint(), SketchLine_579.endPoint())

### Create SketchLine
SketchLine_580 = Sketch_18.addLine(0.6682, -1.50678, 0.6682, -1.50578)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_372).endPoint(), SketchLine_580.startPoint())
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_373).startPoint(), SketchLine_580.endPoint())

### Create SketchLine
SketchLine_581 = Sketch_18.addLine(0.664, -1.55578, 0.664, -1.55728)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_348).endPoint(), SketchLine_581.startPoint())

### Create SketchLine
SketchLine_582 = Sketch_18.addLine(0.664, -1.55728, 0.6682, -1.55728)
Sketch_18.setCoincident(SketchLine_581.endPoint(), SketchLine_582.startPoint())

### Create SketchLine
SketchLine_583 = Sketch_18.addLine(0.6682, -1.55728, 0.6682, -1.55578)
Sketch_18.setCoincident(SketchLine_582.endPoint(), SketchLine_583.startPoint())
Sketch_18.setCoincident(SketchLine_583.endPoint(), SketchArc_323.startPoint())

### Create SketchLine
SketchLine_584 = Sketch_18.addLine(0.6682, -1.50478, 0.6682, -1.5035)
Sketch_18.setCoincident(SketchAPI_Arc(SketchArc_373).endPoint(), SketchLine_584.startPoint())

### Create SketchLine
SketchLine_585 = Sketch_18.addLine(0.6682, -1.5035, 0.664, -1.5035)
Sketch_18.setCoincident(SketchLine_584.endPoint(), SketchLine_585.startPoint())

### Create SketchLine
SketchLine_586 = Sketch_18.addLine(0.664, -1.5035, 0.664, -1.50478)
Sketch_18.setCoincident(SketchLine_585.endPoint(), SketchLine_586.startPoint())
Sketch_18.setCoincident(SketchLine_586.endPoint(), SketchArc_322.startPoint())
model.do()
Sketch_18.setName("BottomStackInsulatingFrame")
Sketch_18.result().setName("BottomStackInsulatingFrame")

### Create Face
Face_18 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "BottomStackInsulatingFrame")])
Face_18.setName("BottomStackInsulatingFrame")
Face_18.result().setName("Face_18_1")

### Create Sketch
Sketch_19 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_587 = Sketch_19.addLine(0.60875, -1.568, 0.60875, -1.5629)

### Create SketchLine
SketchLine_588 = Sketch_19.addLine(0.60875, -1.5629, 0.6155, -1.5629)

### Create SketchLine
SketchLine_589 = Sketch_19.addLine(0.6155, -1.5629, 0.6155, -1.562)

### Create SketchLine
SketchLine_590 = Sketch_19.addLine(0.6155, -1.562, 0.6135, -1.56)

### Create SketchLine
SketchLine_591 = Sketch_19.addLine(0.6135, -1.56, 0.6015, -1.56)

### Create SketchLine
SketchLine_592 = Sketch_19.addLine(0.6015, -1.56, 0.5995, -1.562)

### Create SketchLine
SketchLine_593 = Sketch_19.addLine(0.5995, -1.562, 0.5995, -1.5629)

### Create SketchLine
SketchLine_594 = Sketch_19.addLine(0.5995, -1.5629, 0.6062500000000001, -1.5629)

### Create SketchLine
SketchLine_595 = Sketch_19.addLine(0.6062500000000001, -1.5629, 0.6062500000000001, -1.568)

### Create SketchLine
SketchLine_596 = Sketch_19.addLine(0.6062500000000001, -1.568, 0.60875, -1.568)
Sketch_19.setCoincident(SketchLine_587.endPoint(), SketchLine_588.startPoint())
Sketch_19.setCoincident(SketchLine_588.endPoint(), SketchLine_589.startPoint())
Sketch_19.setCoincident(SketchLine_589.endPoint(), SketchLine_590.startPoint())
Sketch_19.setCoincident(SketchLine_590.endPoint(), SketchLine_591.startPoint())
Sketch_19.setCoincident(SketchLine_591.endPoint(), SketchLine_592.startPoint())
Sketch_19.setCoincident(SketchLine_592.endPoint(), SketchLine_593.startPoint())
Sketch_19.setCoincident(SketchLine_593.endPoint(), SketchLine_594.startPoint())
Sketch_19.setCoincident(SketchLine_594.endPoint(), SketchLine_595.startPoint())
Sketch_19.setCoincident(SketchLine_595.endPoint(), SketchLine_596.startPoint())
Sketch_19.setCoincident(SketchLine_596.endPoint(), SketchLine_587.startPoint())

### Create SketchPoint
SketchPoint_51 = Sketch_19.addPoint(0.60875, -1.56545)
SketchPoint_51.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_51.result(), SketchLine_587.result())
Sketch_19.setFixed(SketchPoint_51.result())

### Create SketchPoint
SketchPoint_52 = Sketch_19.addPoint(0.612125, -1.5629)
SketchPoint_52.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_52.result(), SketchLine_588.result())
Sketch_19.setFixed(SketchPoint_52.result())

### Create SketchPoint
SketchPoint_53 = Sketch_19.addPoint(0.6155, -1.56245)
SketchPoint_53.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_53.result(), SketchLine_589.result())
Sketch_19.setFixed(SketchPoint_53.result())

### Create SketchPoint
SketchPoint_54 = Sketch_19.addPoint(0.6145, -1.561)
SketchPoint_54.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_54.result(), SketchLine_590.result())
Sketch_19.setFixed(SketchPoint_54.result())

### Create SketchPoint
SketchPoint_55 = Sketch_19.addPoint(0.60875, -1.568)

### Create SketchPoint
SketchPoint_56 = Sketch_19.addPoint(0.5277500000000001, -1.568)

### Create SketchMultiTranslation
SketchMultiTranslation_9_objects = [SketchLine_587.result(), SketchLine_588.result(), SketchLine_589.result(), SketchLine_590.result(), SketchLine_591.result(), SketchLine_592.result(), SketchLine_593.result(), SketchLine_594.result(), SketchLine_595.result(), SketchLine_596.result()]
SketchMultiTranslation_9 = Sketch_19.addTranslation(SketchMultiTranslation_9_objects, SketchPoint_55.coordinates(), SketchPoint_56.coordinates(), 8)
[SketchLine_597, SketchLine_598, SketchLine_599, SketchLine_600, SketchLine_601, SketchLine_602, SketchLine_603, SketchLine_604, SketchLine_605, SketchLine_606, SketchLine_607, SketchLine_608, SketchLine_609, SketchLine_610, SketchLine_611, SketchLine_612, SketchLine_613, SketchLine_614, SketchLine_615, SketchLine_616, SketchLine_617, SketchLine_618, SketchLine_619, SketchLine_620, SketchLine_621, SketchLine_622, SketchLine_623, SketchLine_624, SketchLine_625, SketchLine_626, SketchLine_627, SketchLine_628, SketchLine_629, SketchLine_630, SketchLine_631, SketchLine_632, SketchLine_633, SketchLine_634, SketchLine_635, SketchLine_636, SketchLine_637, SketchLine_638, SketchLine_639, SketchLine_640, SketchLine_641, SketchLine_642, SketchLine_643, SketchLine_644, SketchLine_645, SketchLine_646, SketchLine_647, SketchLine_648, SketchLine_649, SketchLine_650, SketchLine_651, SketchLine_652, SketchLine_653, SketchLine_654, SketchLine_655, SketchLine_656, SketchLine_657, SketchLine_658, SketchLine_659, SketchLine_660, SketchLine_661, SketchLine_662, SketchLine_663, SketchLine_664, SketchLine_665, SketchLine_666] = SketchMultiTranslation_9.translatedList()

### Create SketchLine
SketchLine_667 = Sketch_19.addLine(0.6825, -1.56, 0.6805, -1.562)

### Create SketchLine
SketchLine_668 = Sketch_19.addLine(0.6805, -1.562, 0.6805, -1.5629)

### Create SketchLine
SketchLine_669 = Sketch_19.addLine(0.6805, -1.5629, 0.68725, -1.5629)

### Create SketchLine
SketchLine_670 = Sketch_19.addLine(0.68725, -1.5629, 0.68725, -1.568)

### Create SketchLine
SketchLine_671 = Sketch_19.addLine(0.68725, -1.568, 0.6975, -1.568)

### Create SketchLine
SketchLine_672 = Sketch_19.addLine(0.6975, -1.568, 0.6975, -1.56)

### Create SketchLine
SketchLine_673 = Sketch_19.addLine(0.6975, -1.56, 0.69225, -1.56)

### Create SketchLine
SketchLine_674 = Sketch_19.addLine(0.69225, -1.56, 0.69225, -1.563)

### Create SketchLine
SketchLine_675 = Sketch_19.addLine(0.69225, -1.563, 0.69025, -1.563)

### Create SketchLine
SketchLine_676 = Sketch_19.addLine(0.69025, -1.563, 0.69025, -1.56)

### Create SketchLine
SketchLine_677 = Sketch_19.addLine(0.69025, -1.56, 0.6825, -1.56)
Sketch_19.setCoincident(SketchLine_667.endPoint(), SketchLine_668.startPoint())
Sketch_19.setCoincident(SketchLine_668.endPoint(), SketchLine_669.startPoint())
Sketch_19.setCoincident(SketchLine_669.endPoint(), SketchLine_670.startPoint())
Sketch_19.setCoincident(SketchLine_670.endPoint(), SketchLine_671.startPoint())
Sketch_19.setCoincident(SketchLine_671.endPoint(), SketchLine_672.startPoint())
Sketch_19.setCoincident(SketchLine_672.endPoint(), SketchLine_673.startPoint())
Sketch_19.setCoincident(SketchLine_673.endPoint(), SketchLine_674.startPoint())
Sketch_19.setCoincident(SketchLine_674.endPoint(), SketchLine_675.startPoint())
Sketch_19.setCoincident(SketchLine_675.endPoint(), SketchLine_676.startPoint())
Sketch_19.setCoincident(SketchLine_676.endPoint(), SketchLine_677.startPoint())
Sketch_19.setCoincident(SketchLine_677.endPoint(), SketchLine_667.startPoint())

### Create SketchPoint
SketchPoint_57 = Sketch_19.addPoint(0.6815, -1.561)
SketchPoint_57.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_57.result(), SketchLine_667.result())
Sketch_19.setFixed(SketchPoint_57.result())

### Create SketchPoint
SketchPoint_58 = Sketch_19.addPoint(0.6805, -1.56245)
SketchPoint_58.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_58.result(), SketchLine_668.result())
Sketch_19.setFixed(SketchPoint_58.result())

### Create SketchPoint
SketchPoint_59 = Sketch_19.addPoint(0.683875, -1.5629)
SketchPoint_59.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_59.result(), SketchLine_669.result())
Sketch_19.setFixed(SketchPoint_59.result())

### Create SketchPoint
SketchPoint_60 = Sketch_19.addPoint(0.68725, -1.56545)
SketchPoint_60.setAuxiliary(True)
Sketch_19.setCoincident(SketchPoint_60.result(), SketchLine_670.result())
Sketch_19.setFixed(SketchPoint_60.result())
model.do()
Sketch_19.setName("BottomPMTReflectors")
Sketch_19.result().setName("BottomPMTReflectors")

### Create Face
Face_19 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "BottomPMTReflectors")])
Face_19.setName("BottomPMTReflectors")
Face_19.result().setName("Face_19_1")
Face_19.results()[1].setName("Face_19_2")
Face_19.results()[2].setName("Face_19_3")
Face_19.results()[3].setName("Face_19_4")
Face_19.results()[4].setName("Face_19_5")
Face_19.results()[5].setName("Face_19_6")
Face_19.results()[6].setName("Face_19_7")
Face_19.results()[7].setName("Face_19_8")
Face_19.results()[8].setName("Face_19_9")

### Create Sketch
Sketch_20 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_678 = Sketch_20.addLine(0.60875, 0.0819, 0.60875, 0.07679999999999999)

### Create SketchLine
SketchLine_679 = Sketch_20.addLine(0.60875, 0.07679999999999999, 0.6155, 0.07679999999999999)

### Create SketchLine
SketchLine_680 = Sketch_20.addLine(0.6155, 0.07679999999999999, 0.6155, 0.0759)

### Create SketchLine
SketchLine_681 = Sketch_20.addLine(0.6155, 0.0759, 0.6135, 0.07389999999999999)

### Create SketchLine
SketchLine_682 = Sketch_20.addLine(0.6135, 0.07389999999999999, 0.6015, 0.07389999999999999)

### Create SketchLine
SketchLine_683 = Sketch_20.addLine(0.6015, 0.07389999999999999, 0.5995, 0.0759)

### Create SketchLine
SketchLine_684 = Sketch_20.addLine(0.5995, 0.0759, 0.5995, 0.07679999999999999)

### Create SketchLine
SketchLine_685 = Sketch_20.addLine(0.5995, 0.07679999999999999, 0.6062500000000001, 0.07679999999999999)

### Create SketchLine
SketchLine_686 = Sketch_20.addLine(0.6062500000000001, 0.07679999999999999, 0.6062500000000001, 0.0819)

### Create SketchLine
SketchLine_687 = Sketch_20.addLine(0.6062500000000001, 0.0819, 0.60875, 0.0819)
Sketch_20.setCoincident(SketchLine_678.endPoint(), SketchLine_679.startPoint())
Sketch_20.setCoincident(SketchLine_679.endPoint(), SketchLine_680.startPoint())
Sketch_20.setCoincident(SketchLine_680.endPoint(), SketchLine_681.startPoint())
Sketch_20.setCoincident(SketchLine_681.endPoint(), SketchLine_682.startPoint())
Sketch_20.setCoincident(SketchLine_682.endPoint(), SketchLine_683.startPoint())
Sketch_20.setCoincident(SketchLine_683.endPoint(), SketchLine_684.startPoint())
Sketch_20.setCoincident(SketchLine_684.endPoint(), SketchLine_685.startPoint())
Sketch_20.setCoincident(SketchLine_685.endPoint(), SketchLine_686.startPoint())
Sketch_20.setCoincident(SketchLine_686.endPoint(), SketchLine_687.startPoint())
Sketch_20.setCoincident(SketchLine_687.endPoint(), SketchLine_678.startPoint())

### Create SketchPoint
SketchPoint_61 = Sketch_20.addPoint(0.60875, 0.07935)
SketchPoint_61.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_61.result(), SketchLine_678.result())
Sketch_20.setFixed(SketchPoint_61.result())

### Create SketchPoint
SketchPoint_62 = Sketch_20.addPoint(0.612125, 0.07679999999999999)
SketchPoint_62.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_62.result(), SketchLine_679.result())
Sketch_20.setFixed(SketchPoint_62.result())

### Create SketchPoint
SketchPoint_63 = Sketch_20.addPoint(0.6155, 0.07635)
SketchPoint_63.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_63.result(), SketchLine_680.result())
Sketch_20.setFixed(SketchPoint_63.result())

### Create SketchPoint
SketchPoint_64 = Sketch_20.addPoint(0.6145, 0.07489999999999999)
SketchPoint_64.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_64.result(), SketchLine_681.result())
Sketch_20.setFixed(SketchPoint_64.result())

### Create SketchPoint
SketchPoint_65 = Sketch_20.addPoint(0.60875, 0.0819)

### Create SketchPoint
SketchPoint_66 = Sketch_20.addPoint(0.5277500000000001, 0.0819)

### Create SketchMultiTranslation
SketchMultiTranslation_10_objects = [SketchLine_678.result(), SketchLine_679.result(), SketchLine_680.result(), SketchLine_681.result(), SketchLine_682.result(), SketchLine_683.result(), SketchLine_684.result(), SketchLine_685.result(), SketchLine_686.result(), SketchLine_687.result()]
SketchMultiTranslation_10 = Sketch_20.addTranslation(SketchMultiTranslation_10_objects, SketchPoint_65.coordinates(), SketchPoint_66.coordinates(), 8)
[SketchLine_688, SketchLine_689, SketchLine_690, SketchLine_691, SketchLine_692, SketchLine_693, SketchLine_694, SketchLine_695, SketchLine_696, SketchLine_697, SketchLine_698, SketchLine_699, SketchLine_700, SketchLine_701, SketchLine_702, SketchLine_703, SketchLine_704, SketchLine_705, SketchLine_706, SketchLine_707, SketchLine_708, SketchLine_709, SketchLine_710, SketchLine_711, SketchLine_712, SketchLine_713, SketchLine_714, SketchLine_715, SketchLine_716, SketchLine_717, SketchLine_718, SketchLine_719, SketchLine_720, SketchLine_721, SketchLine_722, SketchLine_723, SketchLine_724, SketchLine_725, SketchLine_726, SketchLine_727, SketchLine_728, SketchLine_729, SketchLine_730, SketchLine_731, SketchLine_732, SketchLine_733, SketchLine_734, SketchLine_735, SketchLine_736, SketchLine_737, SketchLine_738, SketchLine_739, SketchLine_740, SketchLine_741, SketchLine_742, SketchLine_743, SketchLine_744, SketchLine_745, SketchLine_746, SketchLine_747, SketchLine_748, SketchLine_749, SketchLine_750, SketchLine_751, SketchLine_752, SketchLine_753, SketchLine_754, SketchLine_755, SketchLine_756, SketchLine_757] = SketchMultiTranslation_10.translatedList()

### Create SketchLine
SketchLine_758 = Sketch_20.addLine(0.6825, 0.07389999999999999, 0.6805, 0.0759)

### Create SketchLine
SketchLine_759 = Sketch_20.addLine(0.6805, 0.0759, 0.6805, 0.07679999999999999)

### Create SketchLine
SketchLine_760 = Sketch_20.addLine(0.6805, 0.07679999999999999, 0.68725, 0.07679999999999999)

### Create SketchLine
SketchLine_761 = Sketch_20.addLine(0.68725, 0.07679999999999999, 0.68725, 0.0819)

### Create SketchLine
SketchLine_762 = Sketch_20.addLine(0.68725, 0.0819, 0.706, 0.0819)

### Create SketchLine
SketchLine_763 = Sketch_20.addLine(0.706, 0.0819, 0.706, 0.07389999999999999)

### Create SketchLine
SketchLine_764 = Sketch_20.addLine(0.706, 0.07389999999999999, 0.6825, 0.07389999999999999)
Sketch_20.setCoincident(SketchLine_758.endPoint(), SketchLine_759.startPoint())
Sketch_20.setCoincident(SketchLine_759.endPoint(), SketchLine_760.startPoint())
Sketch_20.setCoincident(SketchLine_760.endPoint(), SketchLine_761.startPoint())
Sketch_20.setCoincident(SketchLine_761.endPoint(), SketchLine_762.startPoint())
Sketch_20.setCoincident(SketchLine_762.endPoint(), SketchLine_763.startPoint())
Sketch_20.setCoincident(SketchLine_763.endPoint(), SketchLine_764.startPoint())
Sketch_20.setCoincident(SketchLine_764.endPoint(), SketchLine_758.startPoint())

### Create SketchPoint
SketchPoint_67 = Sketch_20.addPoint(0.6815, 0.07489999999999999)
SketchPoint_67.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_67.result(), SketchLine_758.result())
Sketch_20.setFixed(SketchPoint_67.result())

### Create SketchPoint
SketchPoint_68 = Sketch_20.addPoint(0.6805, 0.07635)
SketchPoint_68.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_68.result(), SketchLine_759.result())
Sketch_20.setFixed(SketchPoint_68.result())

### Create SketchPoint
SketchPoint_69 = Sketch_20.addPoint(0.683875, 0.07679999999999999)
SketchPoint_69.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_69.result(), SketchLine_760.result())
Sketch_20.setFixed(SketchPoint_69.result())

### Create SketchPoint
SketchPoint_70 = Sketch_20.addPoint(0.68725, 0.07935)
SketchPoint_70.setAuxiliary(True)
Sketch_20.setCoincident(SketchPoint_70.result(), SketchLine_761.result())
Sketch_20.setFixed(SketchPoint_70.result())
model.do()
Sketch_20.setName("TopPMTReflectors")
Sketch_20.result().setName("TopPMTReflectors")

### Create Face
Face_20 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "TopPMTReflectors")])
Face_20.setName("TopPMTReflectors")
Face_20.result().setName("Face_20_1")
Face_20.results()[1].setName("Face_20_2")
Face_20.results()[2].setName("Face_20_3")
Face_20.results()[3].setName("Face_20_4")
Face_20.results()[4].setName("Face_20_5")
Face_20.results()[5].setName("Face_20_6")
Face_20.results()[6].setName("Face_20_7")
Face_20.results()[7].setName("Face_20_8")
Face_20.results()[8].setName("Face_20_9")

### Create Sketch
Sketch_21 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_765 = Sketch_21.addLine(0.664, -1.5016, 0.664, -0.0007999999999999119)

### Create SketchLine
SketchLine_766 = Sketch_21.addLine(0.664, -0.0007999999999999119, 0.667, -0.0007999999999999119)

### Create SketchLine
SketchLine_767 = Sketch_21.addLine(0.667, -0.0007999999999999119, 0.667, -1.5016)

### Create SketchLine
SketchLine_768 = Sketch_21.addLine(0.667, -1.5016, 0.664, -1.5016)
Sketch_21.setCoincident(SketchLine_765.endPoint(), SketchLine_766.startPoint())
Sketch_21.setCoincident(SketchLine_766.endPoint(), SketchLine_767.startPoint())
Sketch_21.setCoincident(SketchLine_767.endPoint(), SketchLine_768.startPoint())
Sketch_21.setCoincident(SketchLine_768.endPoint(), SketchLine_765.startPoint())

### Create SketchPoint
SketchPoint_71 = Sketch_21.addPoint(0.664, -0.7511999999999999)
SketchPoint_71.setAuxiliary(True)
Sketch_21.setCoincident(SketchPoint_71.result(), SketchLine_765.result())
Sketch_21.setFixed(SketchPoint_71.result())

### Create SketchPoint
SketchPoint_72 = Sketch_21.addPoint(0.6655, -0.0007999999999999119)
SketchPoint_72.setAuxiliary(True)
Sketch_21.setCoincident(SketchPoint_72.result(), SketchLine_766.result())
Sketch_21.setFixed(SketchPoint_72.result())

### Create SketchPoint
SketchPoint_73 = Sketch_21.addPoint(0.667, -0.7511999999999999)
SketchPoint_73.setAuxiliary(True)
Sketch_21.setCoincident(SketchPoint_73.result(), SketchLine_767.result())
Sketch_21.setFixed(SketchPoint_73.result())

### Create SketchPoint
SketchPoint_74 = Sketch_21.addPoint(0.6655, -1.5016)
SketchPoint_74.setAuxiliary(True)
Sketch_21.setCoincident(SketchPoint_74.result(), SketchLine_768.result())
Sketch_21.setFixed(SketchPoint_74.result())
model.do()
Sketch_21.setName("wallPTFE")
Sketch_21.result().setName("wallPTFE")

### Create Face
Face_21 = model.addFace(Part_1_doc, [model.selection("COMPOUND", "wallPTFE")])
Face_21.setName("wallPTFE")
Face_21.result().setName("Face_21_1")

### Create Cut
Cut_1_objects_2 = [model.selection("FACE", "Face_5_1"),
                   model.selection("FACE", "Face_6_1"),
                   model.selection("FACE", "Face_6_2"),
                   model.selection("FACE", "Face_6_3"),
                   model.selection("FACE", "Face_6_4"),
                   model.selection("FACE", "Face_6_5"),
                   model.selection("FACE", "Face_6_6"),
                   model.selection("FACE", "Face_6_7"),
                   model.selection("FACE", "Face_6_8"),
                   model.selection("FACE", "Face_6_9"),
                   model.selection("FACE", "Face_6_10"),
                   model.selection("FACE", "Face_6_11"),
                   model.selection("FACE", "Face_6_12"),
                   model.selection("FACE", "Face_6_13"),
                   model.selection("FACE", "Face_6_14"),
                   model.selection("FACE", "Face_6_15"),
                   model.selection("FACE", "Face_6_16"),
                   model.selection("FACE", "Face_6_17"),
                   model.selection("FACE", "Face_6_18"),
                   model.selection("FACE", "Face_6_19"),
                   model.selection("FACE", "Face_6_20"),
                   model.selection("FACE", "Face_6_21"),
                   model.selection("FACE", "Face_6_22"),
                   model.selection("FACE", "Face_6_23"),
                   model.selection("FACE", "Face_6_24"),
                   model.selection("FACE", "Face_6_25"),
                   model.selection("FACE", "Face_6_26"),
                   model.selection("FACE", "Face_6_27"),
                   model.selection("FACE", "Face_6_28"),
                   model.selection("FACE", "Face_6_29"),
                   model.selection("FACE", "Face_6_30"),
                   model.selection("FACE", "Face_6_31"),
                   model.selection("FACE", "Face_6_32"),
                   model.selection("FACE", "Face_6_33"),
                   model.selection("FACE", "Face_6_34"),
                   model.selection("FACE", "Face_6_35"),
                   model.selection("FACE", "Face_6_36"),
                   model.selection("FACE", "Face_6_37"),
                   model.selection("FACE", "Face_6_38"),
                   model.selection("FACE", "Face_6_39"),
                   model.selection("FACE", "Face_6_40"),
                   model.selection("FACE", "Face_6_41"),
                   model.selection("FACE", "Face_6_42"),
                   model.selection("FACE", "Face_6_43"),
                   model.selection("FACE", "Face_6_44"),
                   model.selection("FACE", "Face_6_45"),
                   model.selection("FACE", "Face_6_46"),
                   model.selection("FACE", "Face_6_47"),
                   model.selection("FACE", "Face_6_48"),
                   model.selection("FACE", "Face_6_49"),
                   model.selection("FACE", "Face_6_50"),
                   model.selection("FACE", "Face_6_51"),
                   model.selection("FACE", "Face_6_52"),
                   model.selection("FACE", "Face_6_53"),
                   model.selection("FACE", "Face_6_54"),
                   model.selection("FACE", "Face_6_55"),
                   model.selection("FACE", "Face_6_56"),
                   model.selection("FACE", "Face_6_57"),
                   model.selection("FACE", "Face_6_58"),
                   model.selection("FACE", "Face_6_59"),
                   model.selection("FACE", "Face_6_60"),
                   model.selection("FACE", "Face_6_61"),
                   model.selection("FACE", "Face_6_62"),
                   model.selection("FACE", "Face_6_63"),
                   model.selection("FACE", "Face_6_64"),
                   model.selection("FACE", "Face_6_65"),
                   model.selection("FACE", "Face_6_66"),
                   model.selection("FACE", "Face_6_67"),
                   model.selection("FACE", "Face_6_68"),
                   model.selection("FACE", "Face_6_69"),
                   model.selection("FACE", "Face_6_70"),
                   model.selection("FACE", "Face_6_71"),
                   model.selection("FACE", "Face_6_72"),
                   model.selection("FACE", "Face_6_73"),
                   model.selection("FACE", "Face_6_74"),
                   model.selection("FACE", "Face_6_75"),
                   model.selection("FACE", "Face_6_76"),
                   model.selection("FACE", "Face_6_77"),
                   model.selection("FACE", "Face_6_78"),
                   model.selection("FACE", "Face_6_79"),
                   model.selection("FACE", "Face_6_80"),
                   model.selection("FACE", "Face_6_81"),
                   model.selection("FACE", "Face_6_82"),
                   model.selection("FACE", "Face_6_83"),
                   model.selection("FACE", "Face_6_84"),
                   model.selection("FACE", "Face_6_85"),
                   model.selection("FACE", "Face_6_86"),
                   model.selection("FACE", "Face_6_87"),
                   model.selection("FACE", "Face_6_88"),
                   model.selection("FACE", "Face_6_89"),
                   model.selection("FACE", "Face_6_90"),
                   model.selection("FACE", "Face_6_91"),
                   model.selection("FACE", "Face_6_92"),
                   model.selection("FACE", "Face_6_93"),
                   model.selection("FACE", "Face_6_94"),
                   model.selection("FACE", "Face_6_95"),
                   model.selection("FACE", "Face_6_96"),
                   model.selection("FACE", "Face_6_97"),
                   model.selection("FACE", "Face_6_98"),
                   model.selection("FACE", "Face_6_99"),
                   model.selection("FACE", "Face_6_100"),
                   model.selection("FACE", "Face_6_101"),
                   model.selection("FACE", "Face_6_102"),
                   model.selection("FACE", "Face_6_103"),
                   model.selection("FACE", "Face_6_104"),
                   model.selection("FACE", "Face_6_105"),
                   model.selection("FACE", "Face_6_106"),
                   model.selection("FACE", "Face_6_107"),
                   model.selection("FACE", "Face_6_108"),
                   model.selection("FACE", "Face_6_109"),
                   model.selection("FACE", "Face_6_110"),
                   model.selection("FACE", "Face_6_111"),
                   model.selection("FACE", "Face_6_112"),
                   model.selection("FACE", "Face_6_113"),
                   model.selection("FACE", "Face_6_114"),
                   model.selection("FACE", "Face_6_115"),
                   model.selection("FACE", "Face_6_116"),
                   model.selection("FACE", "Face_6_117"),
                   model.selection("FACE", "Face_6_118"),
                   model.selection("FACE", "Face_6_119"),
                   model.selection("FACE", "Face_6_120"),
                   model.selection("FACE", "Face_6_121"),
                   model.selection("FACE", "Face_6_122"),
                   model.selection("FACE", "Face_6_123"),
                   model.selection("FACE", "Face_6_124"),
                   model.selection("FACE", "Face_6_125"),
                   model.selection("FACE", "Face_6_126"),
                   model.selection("FACE", "Face_6_127"),
                   model.selection("FACE", "Face_6_128"),
                   model.selection("FACE", "Face_6_129"),
                   model.selection("FACE", "Face_6_130"),
                   model.selection("FACE", "Face_6_131"),
                   model.selection("FACE", "Face_6_132"),
                   model.selection("FACE", "Face_6_133"),
                   model.selection("FACE", "Face_6_134"),
                   model.selection("FACE", "Face_7_1"),
                   model.selection("FACE", "Face_7_2"),
                   model.selection("FACE", "Face_7_3"),
                   model.selection("FACE", "Face_7_4"),
                   model.selection("FACE", "Face_7_5"),
                   model.selection("FACE", "Face_7_6"),
                   model.selection("FACE", "Face_7_7"),
                   model.selection("FACE", "Face_7_8"),
                   model.selection("FACE", "Face_7_9"),
                   model.selection("FACE", "Face_7_10"),
                   model.selection("FACE", "Face_7_11"),
                   model.selection("FACE", "Face_7_12"),
                   model.selection("FACE", "Face_7_13"),
                   model.selection("FACE", "Face_7_14"),
                   model.selection("FACE", "Face_7_15"),
                   model.selection("FACE", "Face_7_16"),
                   model.selection("FACE", "Face_7_17"),
                   model.selection("FACE", "Face_7_18"),
                   model.selection("FACE", "Face_7_19"),
                   model.selection("FACE", "Face_7_20"),
                   model.selection("FACE", "Face_7_21"),
                   model.selection("FACE", "Face_7_22"),
                   model.selection("FACE", "Face_7_23"),
                   model.selection("FACE", "Face_7_24"),
                   model.selection("FACE", "Face_7_25"),
                   model.selection("FACE", "Face_7_26"),
                   model.selection("FACE", "Face_7_27"),
                   model.selection("FACE", "Face_7_28"),
                   model.selection("FACE", "Face_7_29"),
                   model.selection("FACE", "Face_7_30"),
                   model.selection("FACE", "Face_7_31"),
                   model.selection("FACE", "Face_7_32"),
                   model.selection("FACE", "Face_7_33"),
                   model.selection("FACE", "Face_7_34"),
                   model.selection("FACE", "Face_7_35"),
                   model.selection("FACE", "Face_7_36"),
                   model.selection("FACE", "Face_7_37"),
                   model.selection("FACE", "Face_7_38"),
                   model.selection("FACE", "Face_7_39"),
                   model.selection("FACE", "Face_7_40"),
                   model.selection("FACE", "Face_7_41"),
                   model.selection("FACE", "Face_7_42"),
                   model.selection("FACE", "Face_7_43"),
                   model.selection("FACE", "Face_7_44"),
                   model.selection("FACE", "Face_7_45"),
                   model.selection("FACE", "Face_7_46"),
                   model.selection("FACE", "Face_7_47"),
                   model.selection("FACE", "Face_7_48"),
                   model.selection("FACE", "Face_7_49"),
                   model.selection("FACE", "Face_7_50"),
                   model.selection("FACE", "Face_7_51"),
                   model.selection("FACE", "Face_7_52"),
                   model.selection("FACE", "Face_7_53"),
                   model.selection("FACE", "Face_7_54"),
                   model.selection("FACE", "Face_7_55"),
                   model.selection("FACE", "Face_7_56"),
                   model.selection("FACE", "Face_7_57"),
                   model.selection("FACE", "Face_7_58"),
                   model.selection("FACE", "Face_7_59"),
                   model.selection("FACE", "Face_7_60"),
                   model.selection("FACE", "Face_7_61"),
                   model.selection("FACE", "Face_7_62"),
                   model.selection("FACE", "Face_7_63"),
                   model.selection("FACE", "Face_7_64"),
                   model.selection("FACE", "Face_7_65"),
                   model.selection("FACE", "Face_7_66"),
                   model.selection("FACE", "Face_7_67"),
                   model.selection("FACE", "Face_7_68"),
                   model.selection("FACE", "Face_7_69"),
                   model.selection("FACE", "Face_7_70"),
                   model.selection("FACE", "Face_7_71"),
                   model.selection("FACE", "Face_7_72"),
                   model.selection("FACE", "Face_7_73"),
                   model.selection("FACE", "Face_7_74"),
                   model.selection("FACE", "Face_7_75"),
                   model.selection("FACE", "Face_7_76"),
                   model.selection("FACE", "Face_7_77"),
                   model.selection("FACE", "Face_7_78"),
                   model.selection("FACE", "Face_7_79"),
                   model.selection("FACE", "Face_7_80"),
                   model.selection("FACE", "Face_7_81"),
                   model.selection("FACE", "Face_7_82"),
                   model.selection("FACE", "Face_7_83"),
                   model.selection("FACE", "Face_7_84"),
                   model.selection("FACE", "Face_7_85"),
                   model.selection("FACE", "Face_7_86"),
                   model.selection("FACE", "Face_7_87"),
                   model.selection("FACE", "Face_7_88"),
                   model.selection("FACE", "Face_7_89"),
                   model.selection("FACE", "Face_7_90"),
                   model.selection("FACE", "Face_7_91"),
                   model.selection("FACE", "Face_7_92"),
                   model.selection("FACE", "Face_7_93"),
                   model.selection("FACE", "Face_7_94"),
                   model.selection("FACE", "Face_7_95"),
                   model.selection("FACE", "Face_7_96"),
                   model.selection("FACE", "Face_7_97"),
                   model.selection("FACE", "Face_7_98"),
                   model.selection("FACE", "Face_7_99"),
                   model.selection("FACE", "Face_7_100"),
                   model.selection("FACE", "Face_7_101"),
                   model.selection("FACE", "Face_7_102"),
                   model.selection("FACE", "Face_7_103"),
                   model.selection("FACE", "Face_7_104"),
                   model.selection("FACE", "Face_7_105"),
                   model.selection("FACE", "Face_7_106"),
                   model.selection("FACE", "Face_7_107"),
                   model.selection("FACE", "Face_7_108"),
                   model.selection("FACE", "Face_7_109"),
                   model.selection("FACE", "Face_7_110"),
                   model.selection("FACE", "Face_7_111"),
                   model.selection("FACE", "Face_7_112"),
                   model.selection("FACE", "Face_7_113"),
                   model.selection("FACE", "Face_7_114"),
                   model.selection("FACE", "Face_7_115"),
                   model.selection("FACE", "Face_7_116"),
                   model.selection("FACE", "Face_7_117"),
                   model.selection("FACE", "Face_7_118"),
                   model.selection("FACE", "Face_7_119"),
                   model.selection("FACE", "Face_7_120"),
                   model.selection("FACE", "Face_7_121"),
                   model.selection("FACE", "Face_7_122"),
                   model.selection("FACE", "Face_7_123"),
                   model.selection("FACE", "Face_7_124"),
                   model.selection("FACE", "Face_7_125"),
                   model.selection("FACE", "Face_7_126"),
                   model.selection("FACE", "Face_7_127"),
                   model.selection("FACE", "Face_7_128"),
                   model.selection("FACE", "Face_7_129"),
                   model.selection("FACE", "Face_7_130"),
                   model.selection("FACE", "Face_7_131"),
                   model.selection("FACE", "Face_7_132"),
                   model.selection("FACE", "Face_7_133"),
                   model.selection("FACE", "Face_7_134"),
                   model.selection("FACE", "Face_8_1"),
                   model.selection("FACE", "Face_8_2"),
                   model.selection("FACE", "Face_8_3"),
                   model.selection("FACE", "Face_8_4"),
                   model.selection("FACE", "Face_8_5"),
                   model.selection("FACE", "Face_8_6"),
                   model.selection("FACE", "Face_8_7"),
                   model.selection("FACE", "Face_8_8"),
                   model.selection("FACE", "Face_8_9"),
                   model.selection("FACE", "Face_8_10"),
                   model.selection("FACE", "Face_8_11"),
                   model.selection("FACE", "Face_8_12"),
                   model.selection("FACE", "Face_8_13"),
                   model.selection("FACE", "Face_8_14"),
                   model.selection("FACE", "Face_8_15"),
                   model.selection("FACE", "Face_8_16"),
                   model.selection("FACE", "Face_8_17"),
                   model.selection("FACE", "Face_8_18"),
                   model.selection("FACE", "Face_8_19"),
                   model.selection("FACE", "Face_8_20"),
                   model.selection("FACE", "Face_8_21"),
                   model.selection("FACE", "Face_8_22"),
                   model.selection("FACE", "Face_8_23"),
                   model.selection("FACE", "Face_8_24"),
                   model.selection("FACE", "Face_8_25"),
                   model.selection("FACE", "Face_8_26"),
                   model.selection("FACE", "Face_8_27"),
                   model.selection("FACE", "Face_8_28"),
                   model.selection("FACE", "Face_8_29"),
                   model.selection("FACE", "Face_8_30"),
                   model.selection("FACE", "Face_8_31"),
                   model.selection("FACE", "Face_8_32"),
                   model.selection("FACE", "Face_8_33"),
                   model.selection("FACE", "Face_8_34"),
                   model.selection("FACE", "Face_8_35"),
                   model.selection("FACE", "Face_8_36"),
                   model.selection("FACE", "Face_8_37"),
                   model.selection("FACE", "Face_8_38"),
                   model.selection("FACE", "Face_8_39"),
                   model.selection("FACE", "Face_8_40"),
                   model.selection("FACE", "Face_8_41"),
                   model.selection("FACE", "Face_8_42"),
                   model.selection("FACE", "Face_8_43"),
                   model.selection("FACE", "Face_8_44"),
                   model.selection("FACE", "Face_8_45"),
                   model.selection("FACE", "Face_8_46"),
                   model.selection("FACE", "Face_8_47"),
                   model.selection("FACE", "Face_8_48"),
                   model.selection("FACE", "Face_8_49"),
                   model.selection("FACE", "Face_8_50"),
                   model.selection("FACE", "Face_8_51"),
                   model.selection("FACE", "Face_8_52"),
                   model.selection("FACE", "Face_8_53"),
                   model.selection("FACE", "Face_8_54"),
                   model.selection("FACE", "Face_8_55"),
                   model.selection("FACE", "Face_8_56"),
                   model.selection("FACE", "Face_8_57"),
                   model.selection("FACE", "Face_8_58"),
                   model.selection("FACE", "Face_8_59"),
                   model.selection("FACE", "Face_8_60"),
                   model.selection("FACE", "Face_8_61"),
                   model.selection("FACE", "Face_8_62"),
                   model.selection("FACE", "Face_8_63"),
                   model.selection("FACE", "Face_8_64"),
                   model.selection("FACE", "Face_8_65"),
                   model.selection("FACE", "Face_8_66"),
                   model.selection("FACE", "Face_8_67"),
                   model.selection("FACE", "Face_8_68"),
                   model.selection("FACE", "Face_8_69"),
                   model.selection("FACE", "Face_8_70"),
                   model.selection("FACE", "Face_8_71"),
                   model.selection("FACE", "Face_8_72"),
                   model.selection("FACE", "Face_8_73"),
                   model.selection("FACE", "Face_8_74"),
                   model.selection("FACE", "Face_8_75"),
                   model.selection("FACE", "Face_8_76"),
                   model.selection("FACE", "Face_8_77"),
                   model.selection("FACE", "Face_8_78"),
                   model.selection("FACE", "Face_8_79"),
                   model.selection("FACE", "Face_8_80"),
                   model.selection("FACE", "Face_8_81"),
                   model.selection("FACE", "Face_8_82"),
                   model.selection("FACE", "Face_8_83"),
                   model.selection("FACE", "Face_8_84"),
                   model.selection("FACE", "Face_8_85"),
                   model.selection("FACE", "Face_8_86"),
                   model.selection("FACE", "Face_8_87"),
                   model.selection("FACE", "Face_8_88"),
                   model.selection("FACE", "Face_8_89"),
                   model.selection("FACE", "Face_8_90"),
                   model.selection("FACE", "Face_9_1"),
                   model.selection("FACE", "Face_10_1"),
                   model.selection("FACE", "Face_11_1"),
                   model.selection("FACE", "Face_3_1"),
                   model.selection("FACE", "Face_3_2"),
                   model.selection("FACE", "Face_3_3"),
                   model.selection("FACE", "Face_3_4"),
                   model.selection("FACE", "Face_3_5"),
                   model.selection("FACE", "Face_3_6"),
                   model.selection("FACE", "Face_3_7"),
                   model.selection("FACE", "Face_3_8"),
                   model.selection("FACE", "Face_3_9"),
                   model.selection("FACE", "Face_3_10"),
                   model.selection("FACE", "Face_3_11"),
                   model.selection("FACE", "Face_3_12"),
                   model.selection("FACE", "Face_3_13"),
                   model.selection("FACE", "Face_3_14"),
                   model.selection("FACE", "Face_3_15"),
                   model.selection("FACE", "Face_3_16"),
                   model.selection("FACE", "Face_3_17"),
                   model.selection("FACE", "Face_3_18"),
                   model.selection("FACE", "Face_3_19"),
                   model.selection("FACE", "Face_3_20"),
                   model.selection("FACE", "Face_3_21"),
                   model.selection("FACE", "Face_3_22"),
                   model.selection("FACE", "Face_3_23"),
                   model.selection("FACE", "Face_3_24"),
                   model.selection("FACE", "Face_3_25"),
                   model.selection("FACE", "Face_3_26"),
                   model.selection("FACE", "Face_3_27"),
                   model.selection("FACE", "Face_3_28"),
                   model.selection("FACE", "Face_3_29"),
                   model.selection("FACE", "Face_3_30"),
                   model.selection("FACE", "Face_3_31"),
                   model.selection("FACE", "Face_3_32"),
                   model.selection("FACE", "Face_3_33"),
                   model.selection("FACE", "Face_3_34"),
                   model.selection("FACE", "Face_3_35"),
                   model.selection("FACE", "Face_3_36"),
                   model.selection("FACE", "Face_3_37"),
                   model.selection("FACE", "Face_3_38"),
                   model.selection("FACE", "Face_3_39"),
                   model.selection("FACE", "Face_3_40"),
                   model.selection("FACE", "Face_3_41"),
                   model.selection("FACE", "Face_3_42"),
                   model.selection("FACE", "Face_3_43"),
                   model.selection("FACE", "Face_3_44"),
                   model.selection("FACE", "Face_3_45"),
                   model.selection("FACE", "Face_3_46"),
                   model.selection("FACE", "Face_3_47"),
                   model.selection("FACE", "Face_3_48"),
                   model.selection("FACE", "Face_3_49"),
                   model.selection("FACE", "Face_3_50"),
                   model.selection("FACE", "Face_3_51"),
                   model.selection("FACE", "Face_3_52"),
                   model.selection("FACE", "Face_3_53"),
                   model.selection("FACE", "Face_3_54"),
                   model.selection("FACE", "Face_3_55"),
                   model.selection("FACE", "Face_3_56"),
                   model.selection("FACE", "Face_3_57"),
                   model.selection("FACE", "Face_3_58"),
                   model.selection("FACE", "Face_3_59"),
                   model.selection("FACE", "Face_3_60"),
                   model.selection("FACE", "Face_3_61"),
                   model.selection("FACE", "Face_3_62"),
                   model.selection("FACE", "Face_3_63"),
                   model.selection("FACE", "Face_3_64"),
                   model.selection("FACE", "Face_3_65"),
                   model.selection("FACE", "Face_3_66"),
                   model.selection("FACE", "Face_3_67"),
                   model.selection("FACE", "Face_3_68"),
                   model.selection("FACE", "Face_3_69"),
                   model.selection("FACE", "Face_3_70"),
                   model.selection("FACE", "Face_3_71"),
                   model.selection("FACE", "Face_4_1"),
                   model.selection("FACE", "Face_4_2"),
                   model.selection("FACE", "Face_4_3"),
                   model.selection("FACE", "Face_4_4"),
                   model.selection("FACE", "Face_4_5"),
                   model.selection("FACE", "Face_4_6"),
                   model.selection("FACE", "Face_4_7"),
                   model.selection("FACE", "Face_4_8"),
                   model.selection("FACE", "Face_4_9"),
                   model.selection("FACE", "Face_4_10"),
                   model.selection("FACE", "Face_4_11"),
                   model.selection("FACE", "Face_4_12"),
                   model.selection("FACE", "Face_4_13"),
                   model.selection("FACE", "Face_4_14"),
                   model.selection("FACE", "Face_4_15"),
                   model.selection("FACE", "Face_4_16"),
                   model.selection("FACE", "Face_4_17"),
                   model.selection("FACE", "Face_4_18"),
                   model.selection("FACE", "Face_4_19"),
                   model.selection("FACE", "Face_4_20"),
                   model.selection("FACE", "Face_4_21"),
                   model.selection("FACE", "Face_4_22"),
                   model.selection("FACE", "Face_4_23"),
                   model.selection("FACE", "Face_4_24"),
                   model.selection("FACE", "Face_4_25"),
                   model.selection("FACE", "Face_4_26"),
                   model.selection("FACE", "Face_4_27"),
                   model.selection("FACE", "Face_4_28"),
                   model.selection("FACE", "Face_4_29"),
                   model.selection("FACE", "Face_4_30"),
                   model.selection("FACE", "Face_4_31"),
                   model.selection("FACE", "Face_4_32"),
                   model.selection("FACE", "Face_4_33"),
                   model.selection("FACE", "Face_4_34"),
                   model.selection("FACE", "Face_4_35"),
                   model.selection("FACE", "Face_4_36"),
                   model.selection("FACE", "Face_4_37"),
                   model.selection("FACE", "Face_4_38"),
                   model.selection("FACE", "Face_4_39"),
                   model.selection("FACE", "Face_4_40"),
                   model.selection("FACE", "Face_4_41"),
                   model.selection("FACE", "Face_4_42"),
                   model.selection("FACE", "Face_4_43"),
                   model.selection("FACE", "Face_4_44"),
                   model.selection("FACE", "Face_4_45"),
                   model.selection("FACE", "Face_4_46"),
                   model.selection("FACE", "Face_4_47"),
                   model.selection("FACE", "Face_4_48"),
                   model.selection("FACE", "Face_4_49"),
                   model.selection("FACE", "Face_4_50"),
                   model.selection("FACE", "Face_4_51"),
                   model.selection("FACE", "Face_4_52"),
                   model.selection("FACE", "Face_4_53"),
                   model.selection("FACE", "Face_4_54"),
                   model.selection("FACE", "Face_4_55"),
                   model.selection("FACE", "Face_4_56"),
                   model.selection("FACE", "Face_4_57"),
                   model.selection("FACE", "Face_4_58"),
                   model.selection("FACE", "Face_4_59"),
                   model.selection("FACE", "Face_4_60"),
                   model.selection("FACE", "Face_4_61"),
                   model.selection("FACE", "Face_4_62"),
                   model.selection("FACE", "Face_4_63"),
                   model.selection("FACE", "Face_4_64"),
                   model.selection("FACE", "Face_12_1"),
                   model.selection("FACE", "Face_12_2"),
                   model.selection("FACE", "Face_12_3"),
                   model.selection("FACE", "Face_12_4"),
                   model.selection("FACE", "Face_12_5"),
                   model.selection("FACE", "Face_12_6"),
                   model.selection("FACE", "Face_12_7"),
                   model.selection("FACE", "Face_12_8"),
                   model.selection("FACE", "Face_12_9"),
                   model.selection("FACE", "Face_13_1"),
                   model.selection("FACE", "Face_13_2"),
                   model.selection("FACE", "Face_13_3"),
                   model.selection("FACE", "Face_13_4"),
                   model.selection("FACE", "Face_13_5"),
                   model.selection("FACE", "Face_13_6"),
                   model.selection("FACE", "Face_13_7"),
                   model.selection("FACE", "Face_13_8"),
                   model.selection("FACE", "Face_13_9")]
Cut_1 = model.addCut(Part_1_doc, [model.selection("FACE", "Face_2_1"), model.selection("FACE", "Face_1_1")], Cut_1_objects_2, keepSubResults = True)
Cut_1.result().setName("Cut_1_1")
Cut_1.results()[1].setName("Cut_1_2")

### Create Group
Group_1 = model.addGroup(Part_1_doc, "FACE", [model.selection("COMPOUND", "Cut_1_1")])
Group_1.setName("GXeFaces")
Group_1.result().setName("GXeFaces")

### Create Group
Group_2 = model.addGroup(Part_1_doc, "FACE", [model.selection("FACE", "Cut_1_2")])
Group_2.setName("LXeFaces")
Group_2.result().setName("LXeFaces")

### Create GroupShape
GroupShape_1 = model.addGroupShape(Part_1_doc, [model.selection("COMPOUND", "GXeFaces")])
GroupShape_1.setName("GXeGroup")
GroupShape_1.result().setName("GXeGroup")
GroupShape_1.result().subResult(0).setName("GroupShape_1_1_1")
GroupShape_1.result().subResult(1).setName("GroupShape_1_1_2")

### Create GroupShape
GroupShape_2 = model.addGroupShape(Part_1_doc, [model.selection("COMPOUND", "LXeFaces")])
GroupShape_2.setName("LXeGroup")
GroupShape_2.result().setName("LXeGroup")

### Create Group
Group_3_objects = [model.selection("FACE", "Face_14_1"),
                   model.selection("FACE", "Face_15_1"),
                   model.selection("FACE", "Face_16_1"),
                   model.selection("FACE", "Face_17_1"),
                   model.selection("FACE", "Face_18_1"),
                   model.selection("FACE", "Face_19_1"),
                   model.selection("FACE", "Face_19_2"),
                   model.selection("FACE", "Face_19_3"),
                   model.selection("FACE", "Face_19_4"),
                   model.selection("FACE", "Face_19_5"),
                   model.selection("FACE", "Face_19_6"),
                   model.selection("FACE", "Face_19_7"),
                   model.selection("FACE", "Face_19_8"),
                   model.selection("FACE", "Face_19_9"),
                   model.selection("FACE", "Face_20_1"),
                   model.selection("FACE", "Face_20_2"),
                   model.selection("FACE", "Face_20_3"),
                   model.selection("FACE", "Face_20_4"),
                   model.selection("FACE", "Face_20_5"),
                   model.selection("FACE", "Face_20_6"),
                   model.selection("FACE", "Face_20_7"),
                   model.selection("FACE", "Face_20_8"),
                   model.selection("FACE", "Face_20_9"),
                   model.selection("FACE", "Face_21_1")]
Group_3 = model.addGroup(Part_1_doc, "FACE", Group_3_objects)
Group_3.setName("PTFEFaces")
Group_3.result().setName("PTFEFaces")

### Create GroupShape
GroupShape_3 = model.addGroupShape(Part_1_doc, [model.selection("COMPOUND", "PTFEFaces")])
GroupShape_3.setName("PTFE_Group")
GroupShape_3.result().setName("PTFE_Group")
GroupShape_3.result().subResult(0).setName("GroupShape_3_1_1")
GroupShape_3.result().subResult(1).setName("GroupShape_3_1_2")
GroupShape_3.result().subResult(2).setName("GroupShape_3_1_3")
GroupShape_3.result().subResult(3).setName("GroupShape_3_1_4")
GroupShape_3.result().subResult(4).setName("GroupShape_3_1_5")
GroupShape_3.result().subResult(5).setName("GroupShape_3_1_6")
GroupShape_3.result().subResult(6).setName("GroupShape_3_1_7")
GroupShape_3.result().subResult(7).setName("GroupShape_3_1_8")
GroupShape_3.result().subResult(8).setName("GroupShape_3_1_9")
GroupShape_3.result().subResult(9).setName("GroupShape_3_1_10")
GroupShape_3.result().subResult(10).setName("GroupShape_3_1_11")
GroupShape_3.result().subResult(11).setName("GroupShape_3_1_12")
GroupShape_3.result().subResult(12).setName("GroupShape_3_1_13")
GroupShape_3.result().subResult(13).setName("GroupShape_3_1_14")
GroupShape_3.result().subResult(14).setName("GroupShape_3_1_15")
GroupShape_3.result().subResult(15).setName("GroupShape_3_1_16")
GroupShape_3.result().subResult(16).setName("GroupShape_3_1_17")
GroupShape_3.result().subResult(17).setName("GroupShape_3_1_18")
GroupShape_3.result().subResult(18).setName("GroupShape_3_1_19")
GroupShape_3.result().subResult(19).setName("GroupShape_3_1_20")
GroupShape_3.result().subResult(20).setName("GroupShape_3_1_21")
GroupShape_3.result().subResult(21).setName("GroupShape_3_1_22")
GroupShape_3.result().subResult(22).setName("GroupShape_3_1_23")
GroupShape_3.result().subResult(23).setName("GroupShape_3_1_24")

### Create Partition
Partition_1_objects = [model.selection("COMPOUND", "GXeGroup"),
                       model.selection("FACE", "LXeGroup"),
                       model.selection("COMPOUND", "PTFE_Group")]
Partition_1 = model.addPartition(Part_1_doc, Partition_1_objects)
Partition_1.setName("partition_surfaces")
Partition_1.result().setName("partition_surfaces")
Partition_1.result().subResult(0).setName("GXeVol1")
Partition_1.result().subResult(1).setName("TopPMTGuard1")
Partition_1.result().subResult(2).setName("TopPMTGuard2")
Partition_1.result().subResult(3).setName("TopPMTGuard3")
Partition_1.result().subResult(4).setName("TopPMTGuard4")
Partition_1.result().subResult(5).setName("TopPMTGuard5")
Partition_1.result().subResult(6).setName("TopPMTGuard6")
Partition_1.result().subResult(7).setName("TopPMTGuard7")
Partition_1.result().subResult(8).setName("TopPMTGuard8")
Partition_1.result().subResult(9).setName("PTFESpacerGateToAnode")
Partition_1.result().subResult(10).setName("PTFESpacerAnodeToShield")
Partition_1.result().subResult(11).setName("PTFEAboveShield")
Partition_1.result().subResult(12).setName("GXeVol2")
Partition_1.result().subResult(13).setName("GXeVol3")
Partition_1.result().subResult(14).setName("PMTGuard")
Partition_1.result().subResult(15).setName("GXeVol4")
Partition_1.result().subResult(16).setName("LXeVol1")
Partition_1.result().subResult(17).setName("BottomPMTGuard1")
Partition_1.result().subResult(18).setName("BottomPMTGuard2")
Partition_1.result().subResult(19).setName("BottomPMTGuard3")
Partition_1.result().subResult(20).setName("BottomPMTGuard4")
Partition_1.result().subResult(21).setName("BottomPMTGuard5")
Partition_1.result().subResult(22).setName("BottomPMTGuard6")
Partition_1.result().subResult(23).setName("BottomPMTGuard7")
Partition_1.result().subResult(24).setName("BottomPMTGuard8")
Partition_1.result().subResult(25).setName("PTFEBesidesBottomElectrodes")
Partition_1.result().subResult(26).setName("PTFEWall")
Partition_1.result().subResult(27).setName("PTFESplitByInterface")
Partition_1.result().subResult(28).setName("PTFEAboveCathode")
Partition_1.result().subResult(29).setName("LXeVol2")
Partition_1.result().subResult(30).setName("PTFEBelowBottomShield")

### Create Group
Group_4_objects = [model.selection("FACE", "GXeVol1"),
                   model.selection("FACE", "GXeVol2"),
                   model.selection("FACE", "GXeVol3"),
                   model.selection("FACE", "GXeVol4")]
Group_4 = model.addGroup(Part_1_doc, "FACE", Group_4_objects)
Group_4.setName("GXeGroupPostPartition")
Group_4.result().setName("GXeGroupPostPartition")

### Create Group
Group_5 = model.addGroup(Part_1_doc, "FACE", [model.selection("FACE", "LXeVol1"), model.selection("FACE", "LXeVol2")])
Group_5.setName("LXeGroupPostPartition")
Group_5.result().setName("LXeGroupPostPartition")

### Create Group
Group_6_objects = [model.selection("FACE", "PTFESpacerGateToAnode"),
                   model.selection("FACE", "PTFESpacerAnodeToShield"),
                   model.selection("FACE", "PTFEAboveShield"),
                   model.selection("FACE", "PTFEBesidesBottomElectrodes"),
                   model.selection("FACE", "PTFEWall"),
                   model.selection("FACE", "PTFESplitByInterface"),
                   model.selection("FACE", "PTFEAboveCathode"),
                   model.selection("FACE", "PTFEBelowBottomShield"),
                   model.selection("FACE", "TopPMTGuard1"),
                   model.selection("FACE", "TopPMTGuard2"),
                   model.selection("FACE", "TopPMTGuard3"),
                   model.selection("FACE", "TopPMTGuard4"),
                   model.selection("FACE", "TopPMTGuard5"),
                   model.selection("FACE", "TopPMTGuard6"),
                   model.selection("FACE", "TopPMTGuard7"),
                   model.selection("FACE", "TopPMTGuard8"),
                   model.selection("FACE", "BottomPMTGuard1"),
                   model.selection("FACE", "BottomPMTGuard2"),
                   model.selection("FACE", "BottomPMTGuard3"),
                   model.selection("FACE", "BottomPMTGuard4"),
                   model.selection("FACE", "BottomPMTGuard5"),
                   model.selection("FACE", "BottomPMTGuard6"),
                   model.selection("FACE", "BottomPMTGuard7"),
                   model.selection("FACE", "BottomPMTGuard8")]
Group_6 = model.addGroup(Part_1_doc, "FACE", Group_6_objects)
Group_6.setName("PTFE_GroupPostPartition")
Group_6.result().setName("PTFE_GroupPostPartition")

model.end()

###
### SHAPERSTUDY component
###

model.publishToShaperStudy()
import SHAPERSTUDY
Face_14_1, PTFEFaces, = SHAPERSTUDY.shape(model.featureStringId(Face_14))
Face_15_1, PTFEFaces_1, = SHAPERSTUDY.shape(model.featureStringId(Face_15))
Face_16_1, PTFEFaces_2, = SHAPERSTUDY.shape(model.featureStringId(Face_16))
Face_17_1, PTFEFaces_3, = SHAPERSTUDY.shape(model.featureStringId(Face_17))
Face_18_1, PTFEFaces_4, = SHAPERSTUDY.shape(model.featureStringId(Face_18))
Face_19_1, PTFEFaces_5, = SHAPERSTUDY.shape(model.featureStringId(Face_19))
Face_19_2, PTFEFaces_6, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 1))
Face_19_3, PTFEFaces_7, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 2))
Face_19_4, PTFEFaces_8, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 3))
Face_19_5, PTFEFaces_9, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 4))
Face_19_6, PTFEFaces_10, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 5))
Face_19_7, PTFEFaces_11, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 6))
Face_19_8, PTFEFaces_12, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 7))
Face_19_9, PTFEFaces_13, = SHAPERSTUDY.shape(model.featureStringId(Face_19, 8))
Face_20_1, PTFEFaces_14, = SHAPERSTUDY.shape(model.featureStringId(Face_20))
Face_20_2, PTFEFaces_15, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 1))
Face_20_3, PTFEFaces_16, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 2))
Face_20_4, PTFEFaces_17, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 3))
Face_20_5, PTFEFaces_18, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 4))
Face_20_6, PTFEFaces_19, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 5))
Face_20_7, PTFEFaces_20, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 6))
Face_20_8, PTFEFaces_21, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 7))
Face_20_9, PTFEFaces_22, = SHAPERSTUDY.shape(model.featureStringId(Face_20, 8))
Face_21_1, PTFEFaces_23, = SHAPERSTUDY.shape(model.featureStringId(Face_21))
Cut_1_1, GXeFaces, = SHAPERSTUDY.shape(model.featureStringId(Cut_1))
Cut_1_2, LXeFaces, = SHAPERSTUDY.shape(model.featureStringId(Cut_1, 1))
partition_surfaces, GXeGroupPostPartition, LXeGroupPostPartition, PTFE_GroupPostPartition, = SHAPERSTUDY.shape(model.featureStringId(Partition_1))
###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(partition_surfaces)
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
GXeGroupPostPartition_1 = Mesh_1.GroupOnGeom(GXeGroupPostPartition,'GXeGroupPostPartition',SMESH.FACE)
LXeGroupPostPartition_1 = Mesh_1.GroupOnGeom(LXeGroupPostPartition,'LXeGroupPostPartition',SMESH.FACE)
PTFE_GroupPostPartition_1 = Mesh_1.GroupOnGeom(PTFE_GroupPostPartition,'PTFE_GroupPostPartition',SMESH.FACE)
isDone = Mesh_1.Compute()
[ GXeGroupPostPartition_1, LXeGroupPostPartition_1, PTFE_GroupPostPartition_1 ] = Mesh_1.GetGroups()


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(GXeGroupPostPartition_1, 'GXeGroupPostPartition')
smesh.SetName(PTFE_GroupPostPartition_1, 'PTFE_GroupPostPartition')
smesh.SetName(LXeGroupPostPartition_1, 'LXeGroupPostPartition')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
