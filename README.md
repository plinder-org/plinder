# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/plinder-org/plinder/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                                           |    Stmts |     Miss |   Branch |   BrPart |        Cover |   Missing |
|--------------------------------------------------------------- | -------: | -------: | -------: | -------: | -----------: | --------: |
| src//plinder/\_\_init\_\_.py                                   |        4 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/\_version.py                                      |       15 |        4 |        0 |        0 |     73.3333% |     16-19 |
| src//plinder/core/\_\_init\_\_.py                              |        5 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/index/\_\_init\_\_.py                        |        2 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/index/system.py                              |      197 |       55 |       96 |       27 |     64.5051% |13, 40, 59->58, 73-74, 78->77, 96->95, 110, 114->113, 127->126, 140->139, 153->152, 166->165, 179->178, 195->194, 209-210, 214->213, 230->229, 243->242, 255-258, 261->260, 270-279, 297-321, 324->323, 335-336, 340->339, 350-351, 361->360, 372->371, 383->382, 389-403, 406->405, 420->419, 422, 427->432, 428->427, 433 |
| src//plinder/core/index/utils.py                               |      133 |       18 |       68 |       13 |     82.5871% |30->29, 75, 112->111, 162, 170-175, 177-178, 180-181, 199, 252, 260, 275, 291, 293 |
| src//plinder/core/loader/\_\_init\_\_.py                       |        2 |        2 |        0 |        0 |      0.0000% |       4-6 |
| src//plinder/core/loader/dataset.py                            |       31 |       31 |        4 |        0 |      0.0000% |      3-91 |
| src//plinder/core/loader/featurizer.py                         |       29 |       29 |        0 |        0 |      0.0000% |     1-102 |
| src//plinder/core/loader/transforms.py                         |       13 |       13 |        0 |        0 |      0.0000% |      1-20 |
| src//plinder/core/loader/utils.py                              |       56 |       56 |       12 |        0 |      0.0000% |     1-210 |
| src//plinder/core/scores/\_\_init\_\_.py                       |        8 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/scores/clusters.py                           |       19 |        2 |        4 |        2 |     82.6087% |19->18, 51-52 |
| src//plinder/core/scores/index.py                              |       27 |        4 |       10 |        4 |     78.3784% |41, 43, 53, 58 |
| src//plinder/core/scores/ligand.py                             |       46 |        3 |       14 |        5 |     86.6667% |19->18, 49-52, 58->57, 66, 93->92 |
| src//plinder/core/scores/links.py                              |       32 |        7 |       16 |        6 |     68.7500% |20->19, 42->exit, 44-45, 47-48, 51, 65-66 |
| src//plinder/core/scores/protein.py                            |       77 |        6 |       38 |       12 |     84.3478% |20->19, 47->53, 61-62, 67->66, 94->93, 126->125, 160, 164->174, 170, 173, 181->183, 190 |
| src//plinder/core/scores/query.py                              |      104 |        9 |       56 |        8 |     88.1250% |49, 97, 104, 173->181, 176-178, 216, 218, 264 |
| src//plinder/core/split/\_\_init\_\_.py                        |        2 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/split/plot.py                                |      302 |      201 |      104 |       11 |     33.4975% |19-20, 31-39, 45-57, 142->141, 156, 158, 159->161, 171, 176, 189-205, 258-261, 295-313, 329-330, 364, 369-397, 400-428, 431-473, 476-603, 606-614, 617-622, 625-626, 629-660, 663-710, 713-755, 762, 764-766, 768-773, 777-820, 832 |
| src//plinder/core/split/utils.py                               |       19 |        0 |        4 |        1 |     95.6522% |    17->16 |
| src//plinder/core/structure/\_\_init\_\_.py                    |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/structure/atoms.py                           |      175 |       83 |       60 |        7 |     46.3830% |85, 99-101, 128-132, 134, 153-156, 160, 211-241, 256-279, 285-289, 293-295, 299, 305-306, 312, 320-337, 345-346, 355-356, 364-375, 382-389 |
| src//plinder/core/structure/contacts.py                        |       88 |        8 |       30 |        5 |     88.9831% |93->95, 112-113, 191->195, 207, 220, 226-232, 258-261 |
| src//plinder/core/structure/diffdock\_utils.py                 |      177 |      146 |       60 |        0 |     13.0802% |26, 32, 38-43, 58-79, 86-93, 96-98, 102-121, 127-169, 173-188, 372-375, 379-446, 450-464, 476-520 |
| src//plinder/core/structure/models.py                          |       12 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/structure/structure.py                       |      327 |      182 |      140 |       39 |     40.6852% |37, 63-73, 82, 86-89, 158->157, 162, 172, 178-181, 198-206, 210, 212, 214, 237-238, 241-262, 279-281, 290-310, 316-319, 329-390, 393-394, 406, 408, 445->444, 449-455, 458->457, 460-463, 466->465, 468-473, 476->475, 480-491, 494->493, 498-504, 514->513, 519, 522->521, 524-533, 536->535, 538-543, 546->545, 548-553, 556->555, 562->561, 563-573, 576->575, 578-579, 582->581, 584-590, 593->592, 595-596, 599->598, 610->609, 612-618, 621->620, 623-625, 628->627, 630-632, 635->634, 637-639, 642->641, 644-650, 653->652, 655-658, 661->660, 663-669, 672->671, 674-678, 681->680, 683-687, 690->689, 692-696, 699->698, 701-705, 708->707, 711-716, 719->718, 720 |
| src//plinder/core/structure/superimpose.py                     |       57 |        7 |       16 |        5 |     80.8219% |92, 104-107, 115, 138, 162->166 |
| src//plinder/core/structure/surgery.py                         |       65 |       51 |       18 |        0 |     16.8675% |19-30, 38-48, 54-86, 97, 118-129 |
| src//plinder/core/structure/vendored.py                        |      367 |      216 |      134 |       20 |     37.3253% |30-32, 36-38, 42->41, 44-53, 70-72, 77-87, 96, 103-104, 124-127, 137->139, 140, 141->143, 161-162, 167-173, 179-184, 194-284, 340->exit, 342, 352, 377->379, 379->381, 384, 424->427, 427->433, 440->442, 442->444, 444->439, 470-472, 483-503, 515, 517, 531-551, 558-578, 613-638, 671-688, 694-707, 763-798, 817-826, 842-844, 862-865, 885-903 |
| src//plinder/core/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/utils/config.py                              |      125 |        4 |       44 |        5 |     94.6746% |130->129, 134, 155, 262, 284 |
| src//plinder/core/utils/constants.py                           |       48 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/utils/cpl.py                                 |      109 |       37 |       46 |       12 |     59.3548% |28->27, 29-37, 45->44, 46, 50->49, 51, 58, 69, 73->72, 74-85, 89->88, 102-105, 114, 152->154, 159-168, 178 |
| src//plinder/core/utils/dataclass.py                           |       33 |       26 |       24 |        0 |     12.2807% |10-18, 29-74 |
| src//plinder/core/utils/dec.py                                 |       19 |        0 |        2 |        1 |     95.2381% |    18->17 |
| src//plinder/core/utils/gcs.py                                 |       93 |       25 |       40 |       13 |     65.4135% |26->25, 28->50, 33-35, 47, 54-58, 66->65, 79->78, 93->92, 109->108, 117-120, 124->123, 132-133, 135, 145, 157-164 |
| src//plinder/core/utils/io.py                                  |       47 |        7 |       14 |        5 |     80.3279% |27->26, 47->46, 54->62, 57-58, 66->65, 67-81 |
| src//plinder/core/utils/load\_systems.py                       |       23 |       23 |       12 |        0 |      0.0000% |      4-41 |
| src//plinder/core/utils/log.py                                 |       37 |       12 |       14 |        2 |     64.7059% |12-13, 63-66, 83-88 |
| src//plinder/core/utils/schemas.py                             |        9 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/utils/unpack.py                              |       88 |        8 |       48 |        7 |     88.9706% |33, 39, 45, 114, 127->132, 134->133, 138-140, 150 |
| src//plinder/data/\_\_init\_\_.py                              |        6 |        2 |        0 |        0 |     66.6667% |       8-9 |
| src//plinder/data/\_version.py                                 |       15 |       15 |        0 |        0 |      0.0000% |     18-40 |
| src//plinder/data/clusters.py                                  |      118 |       37 |       26 |        6 |     65.9722% |195-232, 257-258, 267->exit, 274, 305-306, 320-321, 332-369, 380 |
| src//plinder/data/column\_descriptions/\_\_init\_\_.py         |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/common/\_\_init\_\_.py                       |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/common/\_version.py                          |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/common/constants.py                          |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/common/log.py                                |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/databases.py                                 |       59 |        3 |       26 |        2 |     94.1176% |13, 137->148, 153-154 |
| src//plinder/data/docs.py                                      |       57 |        0 |       24 |        1 |     98.7654% |    85->87 |
| src//plinder/data/final\_structure\_qc.py                      |      152 |       34 |       30 |        6 |     76.9231% |22, 45-47, 69-72, 114-123, 144-145, 147, 151-152, 154, 252-253, 280-283, 367-372 |
| src//plinder/data/get\_system\_annotations.py                  |      106 |       46 |       44 |        5 |     47.3333% |59, 62-63, 73-112, 141->145, 154->exit, 158-175, 187-188 |
| src//plinder/data/leakage.py                                   |       93 |       83 |       16 |        0 |      9.1743% |14-15, 27-51, 59-90, 100-163, 177-229 |
| src//plinder/data/pipeline/\_\_init\_\_.py                     |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/pipeline/config.py                           |      117 |        3 |       42 |        3 |     96.2264% |160, 183, 187 |
| src//plinder/data/pipeline/io.py                               |      272 |       60 |      126 |       22 |     71.8593% |38->37, 77->76, 119-145, 153->152, 197->196, 236->235, 269->277, 281->280, 316->315, 363->361, 413->412, 433, 439, 444, 448, 462->461, 469, 474-476, 500-505, 515-524, 528->527, 547-576, 580->579 |
| src//plinder/data/pipeline/mpqueue.py                          |       65 |       15 |       18 |        6 |     74.6988% |27, 32-35, 51-56, 63->65, 79->78, 80-85, 88->exit, 135->134 |
| src//plinder/data/pipeline/pipeline.py                         |      181 |       39 |       92 |       40 |     70.3297% |53-57, 60->59, 69->68, 76->75, 84->83, 91->90, 105->104, 120->119, 127->126, 136->135, 143->142, 153->152, 163->162, 167->166, 174->173, 183->182, 193->192, 194, 197->196, 206->205, 213->212, 220->219, 230->229, 238->237, 246->245, 247-250, 258->257, 259-260, 263->262, 264, 267->266, 268-280, 283->282, 286-289, 296->295, 297-301, 304->303, 305, 308->307, 309-313, 316->315, 317, 323->322, 324-330, 333->332, 334, 340->339, 341-345, 348->347, 349, 355->354, 356-360, 368->367, 369-374, 377->376, 378-389, 392->391, 393-394, 440-441 |
| src//plinder/data/pipeline/tasks.py                            |      379 |      151 |      126 |       12 |     56.2376% |161, 174->177, 178-179, 218->227, 287-288, 297-298, 327-328, 334, 341, 372, 403->405, 409-413, 414->416, 424-427, 428->397, 466, 723-731, 737-738, 753-772, 793-812, 821-822, 844-864, 874-914, 922-923, 931-947, 955-977, 985-986, 999-1000, 1013, 1022-1028, 1044-1070, 1079-1090, 1101-1125 |
| src//plinder/data/pipeline/transform.py                        |      140 |       81 |       85 |        4 |     39.1111% |26, 29-30, 34, 48, 125-169, 187-210, 234-319, 356 |
| src//plinder/data/pipeline/utils.py                            |      375 |      134 |      154 |       17 |     61.6257% |28, 41->40, 50-53, 77->76, 103-104, 110->109, 142-143, 155-156, 165, 167->170, 171-172, 174-175, 206-207, 229->227, 293, 324-331, 421->420, 479-500, 509, 568->565, 590, 594->596, 605-610, 621-651, 662-686, 705-736, 751-752, 759-767, 780-797, 811-820 |
| src//plinder/data/save\_linked\_structures.py                  |      163 |       99 |       64 |        5 |     35.6828% |23-32, 37-44, 70-108, 136-137, 148, 180->182, 187-188, 223-224, 236-246, 258-290, 301-327, 340-347, 383-385 |
| src//plinder/data/smallmolecules.py                            |      110 |       36 |       14 |        6 |     64.5161% |19-20, 26-27, 33-34, 40-41, 47-48, 54-55, 61-62, 71-72, 83-89, 101, 105-106, 116-118, 125->124, 127-129, 139-140, 152->154, 163-166 |
| src//plinder/data/splits.py                                    |      296 |      231 |       91 |        3 |     18.6047% |33->exit, 90->exit, 97->exit, 121-126, 144-151, 155-162, 166-169, 191-199, 210-222, 242-331, 352-392, 416-464, 497-515, 557-615, 652-897, 901-961 |
| src//plinder/data/structure/\_\_init\_\_.py                    |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/structure/atoms.py                           |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/structure/contacts.py                        |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/utils/annotations/\_\_init\_\_.py            |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/utils/annotations/aggregate\_annotations.py  |      593 |       79 |      387 |       71 |     80.8163% |125->124, 138->137, 151->150, 158->157, 165->164, 172->171, 179->178, 186->185, 193->192, 200->199, 207->206, 214->213, 228->227, 234->exit, 237, 240->239, 247->246, 254->253, 266->265, 270-281, 284->283, 288-295, 307, 315, 323->322, 330->329, 343->342, 350->349, 452->451, 470->469, 480->479, 487->486, 492, 500->499, 505, 513->512, 522->521, 533->532, 539, 547->546, 552, 560->559, 565, 593->591, 632-633, 664, 670, 678-682, 691-692, 710, 721, 733, 833->835, 846->845, 868->867, 934->939, 937-938, 971->973, 975->979, 994->996, 1020->1032, 1067->1056, 1087->1086, 1091, 1124-1138, 1141, 1146, 1149, 1184, 1276-1277, 1339-1342, 1345-1346, 1352->exit, 1395, 1418, 1422, 1434-1445, 1454-1475 |
| src//plinder/data/utils/annotations/get\_ligand\_validation.py |      137 |       12 |       30 |       10 |     86.8263% |59->58, 66-70, 82, 123-127, 187->186, 196, 299->298, 307-308, 311, 330->332, 331->330, 332->331, 338 |
| src//plinder/data/utils/annotations/get\_similarity\_scores.py |      519 |      408 |      268 |       10 |     17.9161% |81-129, 133-136, 150, 207->220, 221, 227->230, 232, 233->230, 242-279, 312, 315->314, 322-323, 365-375, 383-391, 396-403, 412-497, 504-528, 533-596, 599-640, 643-650, 664-744, 755-809, 814-819, 824-900, 903-906, 909-918, 925-986, 991-1039 |
| src//plinder/data/utils/annotations/interaction\_utils.py      |      148 |        6 |       82 |        5 |     93.4783% |109-110, 158-161, 213, 322, 351->359 |
| src//plinder/data/utils/annotations/interface\_gap.py          |       87 |        8 |       28 |        6 |     87.8261% |53, 119, 128, 130, 166-167, 185-186 |
| src//plinder/data/utils/annotations/ligand\_utils.py           |      588 |       59 |      278 |       51 |     86.1432% |120->119, 143->exit, 175, 184, 199, 245-247, 285-290, 312->301, 320, 342, 344, 346, 349-352, 356->355, 443->442, 458->457, 468->467, 516-518, 562, 567, 593-594, 606-607, 618-619, 630-631, 643-644, 655-656, 886->891, 893-896, 928-931, 959->958, 1017, 1141-1142, 1144, 1164->1163, 1170->1172, 1175->1174, 1186->1185, 1195->1194, 1204->1203, 1211->1210, 1221->1220, 1231->1230, 1237->1239, 1243, 1265->1264, 1272->1271, 1282->1281, 1287, 1291->1290, 1298->1297, 1305->1304, 1312->1311, 1316-1323, 1326->1325, 1332-1333, 1337->1336, 1344-1345, 1378, 1435, 1472->1474, 1527->1531 |
| src//plinder/data/utils/annotations/mmpdb\_utils.py            |      124 |       66 |       40 |        0 |     45.1220% |35, 47-165, 173-185, 405-414 |
| src//plinder/data/utils/annotations/protein\_utils.py          |      138 |        4 |       68 |       12 |     92.2330% |61-62, 88->86, 95->102, 97->95, 99->95, 153->155, 262->261, 266, 304->303, 337->339, 350->349, 357->356, 361, 388->387 |
| src//plinder/data/utils/annotations/rdkit\_utils.py            |      265 |       78 |       78 |       13 |     68.2216% |53-61, 103, 105, 110-115, 126->129, 131->118, 159-165, 169-176, 204-205, 288-289, 313->346, 322->337, 338-343, 349, 362-369, 391-392, 442-443, 457-466, 471, 473-477, 488, 491->500, 496-497, 505-532 |
| src//plinder/data/utils/annotations/save\_utils.py             |       80 |        4 |       38 |        2 |     93.2203% |88-90, 148 |
| src//plinder/data/utils/annotations/utils.py                   |       42 |        2 |       30 |        6 |     88.8889% |14->13, 31->33, 37->36, 51, 53, 57->56 |
| src//plinder/data/utils/cluster.py                             |       43 |       43 |       12 |        0 |      0.0000% |     3-101 |
| src//plinder/data/utils/diffdock\_utils.py                     |      177 |      177 |       60 |        0 |      0.0000% |     6-520 |
| src//plinder/data/utils/tanimoto.py                            |      105 |       11 |       26 |        8 |     85.4962% |19, 31-32, 56-57, 62, 75, 137, 192, 194, 196 |
| src//plinder/eval/\_\_init\_\_.py                              |        5 |        2 |        0 |        0 |     60.0000% |       7-8 |
| src//plinder/eval/docking/\_\_init\_\_.py                      |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/eval/docking/make\_plots.py                       |      126 |        5 |       30 |        6 |     92.9487% |43->42, 75->78, 84-85, 119, 264, 275 |
| src//plinder/eval/docking/stratify\_test\_set.py               |      109 |       11 |       42 |        9 |     85.4305% |64->71, 126->125, 135, 190->189, 233-236, 238-241, 274-277, 323, 335 |
| src//plinder/eval/docking/utils.py                             |      203 |       23 |       78 |       16 |     83.9858% |57->56, 76->75, 80, 86, 123->122, 166->168, 172->171, 211-227, 326->315, 333->365, 357-362, 365->370, 378, 383->381, 387, 398, 410-419, 441->443, 443->458 |
| src//plinder/eval/docking/write\_scores.py                     |      101 |       49 |       45 |        3 |     45.8904% |56-86, 126-171, 213->218, 242-243, 295, 313 |
| src//plinder/methods/\_\_init\_\_.py                           |        0 |        0 |        0 |        0 |    100.0000% |           |
|                                                      **TOTAL** | **9328** | **3403** | **3726** |  **576** | **60.1884%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/plinder-org/plinder/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/plinder-org/plinder/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/plinder-org/plinder/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/plinder-org/plinder/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Fplinder-org%2Fplinder%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/plinder-org/plinder/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.