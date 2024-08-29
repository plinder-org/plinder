# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/plinder-org/plinder/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                                          |    Stmts |     Miss |   Branch |   BrPart |        Cover |   Missing |
|-------------------------------------------------------------- | -------: | -------: | -------: | -------: | -----------: | --------: |
| src/plinder/\_\_init\_\_.py                                   |        4 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/\_version.py                                      |       15 |        4 |        0 |        0 |     73.3333% |     16-19 |
| src/plinder/core/\_\_init\_\_.py                              |       13 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/index/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/index/system.py                              |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/index/utils.py                               |       64 |        5 |       24 |        8 |     85.2273% |25->24, 76, 81-83, 120->119, 158->157, 163->158, 164->163, 170 |
| src/plinder/core/loader/\_\_init\_\_.py                       |        2 |        1 |        0 |        0 |     50.0000% |         6 |
| src/plinder/core/loader/loader.py                             |       37 |       34 |       18 |        0 |      5.4545% |      8-92 |
| src/plinder/core/scores/\_\_init\_\_.py                       |        8 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/scores/clusters.py                           |       19 |        2 |        4 |        2 |     82.6087% |19->18, 51-52 |
| src/plinder/core/scores/index.py                              |       16 |        1 |        2 |        1 |     88.8889% |        39 |
| src/plinder/core/scores/ligand.py                             |       46 |        3 |       14 |        5 |     86.6667% |19->18, 49-52, 58->57, 66, 93->92 |
| src/plinder/core/scores/links.py                              |       19 |        0 |        4 |        1 |     95.6522% |    20->19 |
| src/plinder/core/scores/protein.py                            |       50 |        4 |       22 |        6 |     86.1111% |19->18, 46->52, 60-61, 66->65, 76-77, 93->92 |
| src/plinder/core/scores/query.py                              |      104 |        9 |       56 |        8 |     88.1250% |49, 73, 97, 104, 173->181, 176-178, 216, 264 |
| src/plinder/core/split/\_\_init\_\_.py                        |        2 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/split/plot.py                                |      302 |      190 |      104 |       12 |     36.4532% |17-18, 38, 44-56, 141->140, 155, 157, 158->160, 170, 175, 188-204, 257-260, 294-312, 328-329, 378-396, 399-430, 433-475, 478-605, 608-616, 627-628, 631-662, 665-712, 715-757, 764, 766-768, 771-775, 779-822, 834 |
| src/plinder/core/split/utils.py                               |       59 |       10 |       14 |        3 |     79.4521% |19->18, 37-44, 103->105, 181-182 |
| src/plinder/core/system/\_\_init\_\_.py                       |        2 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/system/system.py                             |      119 |       23 |       58 |       18 |     70.0565% |29, 48->47, 62-63, 67->66, 82, 85->84, 99->102, 100->102, 105->104, 118->117, 131->130, 144->143, 157->156, 168->170, 173->172, 184->186, 189->188, 205->204, 218->217, 230-247, 250->249, 259-260, 278-289 |
| src/plinder/core/system/utils.py                              |       23 |       23 |       12 |        0 |      0.0000% |      3-39 |
| src/plinder/core/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/utils/config.py                              |      123 |        4 |       44 |        5 |     94.6108% |125->124, 129, 150, 259, 281 |
| src/plinder/core/utils/cpl.py                                 |      102 |       43 |       42 |       11 |     48.6111% |26->25, 27-35, 43->42, 44, 48->47, 49, 56, 61-67, 71->70, 72-83, 87->86, 100-103, 112, 139->141, 145-154, 165 |
| src/plinder/core/utils/dec.py                                 |       19 |        0 |        2 |        1 |     95.2381% |    18->17 |
| src/plinder/core/utils/gcs.py                                 |       93 |       25 |       40 |       14 |     64.6617% |26->25, 28->50, 33-35, 47, 54-58, 66->65, 79->78, 93->92, 109->108, 117-120, 124->123, 132-133, 135, 142->exit, 145, 157-164 |
| src/plinder/core/utils/log.py                                 |       37 |       12 |       14 |        2 |     64.7059% |12-13, 63-66, 83-88 |
| src/plinder/core/utils/schemas.py                             |        9 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/utils/unpack.py                              |       54 |        3 |       28 |        4 |     91.4634% |28, 34, 40, 80->85 |
| src/plinder/data/\_\_init\_\_.py                              |        9 |        4 |        0 |        0 |     55.5556% |7-8, 20-21 |
| src/plinder/data/\_version.py                                 |       15 |       15 |        0 |        0 |      0.0000% |     18-40 |
| src/plinder/data/clusters.py                                  |      118 |       37 |       26 |        6 |     65.9722% |194-231, 256-257, 266->exit, 273, 304-305, 319-320, 331-368, 379 |
| src/plinder/data/common/\_\_init\_\_.py                       |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/common/\_version.py                          |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/common/constants.py                          |       30 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/common/log.py                                |        3 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/databases.py                                 |       59 |        3 |       26 |        6 |     89.4118% |13, 95->103, 96->95, 97->96, 137->148, 153-154, 221->exit |
| src/plinder/data/final\_structure\_qc.py                      |      144 |       30 |       28 |        7 |     77.3256% |22, 45-47, 69-72, 114-123, 147, 245-246, 273-276, 360-365, 457, 480->482 |
| src/plinder/data/get\_system\_annotations.py                  |      106 |       46 |       44 |        6 |     46.6667% |59, 62-63, 73-112, 141->145, 154->exit, 156->exit, 158-175, 187-188 |
| src/plinder/data/leakage.py                                   |       93 |       83 |       16 |        0 |      9.1743% |14-15, 27-51, 59-90, 100-163, 177-229 |
| src/plinder/data/pipeline/\_\_init\_\_.py                     |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/pipeline/config.py                           |      109 |        4 |       40 |        3 |     95.3020% |151, 174, 178, 215 |
| src/plinder/data/pipeline/io.py                               |      303 |       62 |      138 |       34 |     71.4286% |40->39, 60->59, 90->95, 93->95, 99->98, 141-167, 169->171, 175->174, 207->210, 208->207, 219->218, 258->257, 291->299, 303->302, 332->334, 338->337, 385->383, 435->434, 455, 461, 466, 470, 484->483, 491, 496-498, 522-527, 537-546, 550->549, 557->565, 560-561, 563->565, 569->568, 588-617, 621->620, 649->651 |
| src/plinder/data/pipeline/mpqueue.py                          |       65 |       15 |       18 |        8 |     72.2892% |27, 32-35, 51-56, 63->65, 68->exit, 79->78, 80-85, 88->exit, 132->exit, 135->134 |
| src/plinder/data/pipeline/pipeline.py                         |      160 |       28 |       82 |       35 |     73.1405% |53-57, 60->59, 69->68, 76->75, 84->83, 91->90, 105->104, 120->119, 127->126, 136->135, 143->142, 153->152, 163->162, 167->166, 174->173, 183->182, 193->192, 194, 197->196, 206->205, 213->212, 220->219, 230->229, 238->237, 246->245, 247-250, 258->257, 259-260, 263->262, 264, 267->266, 268-280, 283->282, 286-289, 296->295, 297-301, 304->303, 305, 308->307, 309-313, 316->315, 317, 323->322, 324-330, 333->332, 334, 340->339, 341, 392-393 |
| src/plinder/data/pipeline/tasks.py                            |      334 |      110 |      104 |       22 |     61.6438% |145->exit, 156->145, 159, 172->175, 176-177, 216->225, 279->312, 280->279, 285-286, 295-296, 325-326, 332, 337->340, 338->337, 339, 370, 394->393, 395->394, 401->403, 407-411, 412->414, 422-425, 426->395, 464, 483->479, 484->483, 721-729, 735-736, 751-770, 791-810, 819-820, 842-862, 872-912, 920-921, 929-945, 953-975, 983-984, 997-998, 1012-1029 |
| src/plinder/data/pipeline/transform.py                        |      140 |       81 |       85 |        5 |     38.6667% |26, 29-30, 34, 48, 99->101, 125-169, 187-210, 234-319, 356 |
| src/plinder/data/pipeline/utils.py                            |      301 |       83 |      126 |       21 |     69.0867% |26-27, 40->39, 49-52, 76->75, 102-103, 109->108, 141-142, 143->139, 147->143, 149->147, 154-155, 164, 166->169, 170-171, 173-174, 205-206, 228->226, 292, 323-330, 420->419, 467->465, 473, 475->491, 485-488, 491->510, 504-507, 511, 521-540, 557-586, 601-602, 617-629, 641-650 |
| src/plinder/data/save\_linked\_structures.py                  |      154 |      154 |       62 |        0 |      0.0000% |     3-362 |
| src/plinder/data/smallmolecules.py                            |      110 |       41 |       14 |        5 |     59.6774% |19-20, 26-27, 33-34, 40-41, 47-48, 54-55, 61-62, 71-72, 83-89, 101, 105-106, 116-118, 125->124, 127-129, 139-140, 152-156, 163-166 |
| src/plinder/data/splits.py                                    |      298 |      232 |       91 |        3 |     18.7661% |34->exit, 91->exit, 98->exit, 122-127, 145-152, 156-163, 167-170, 192-200, 211-223, 243-333, 354-394, 418-466, 499-517, 559-617, 654-899, 903-963 |
| src/plinder/data/structure/\_\_init\_\_.py                    |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/structure/atoms.py                           |      301 |      148 |       78 |       20 |     47.2296% |73-75, 85-96, 103, 109->113, 115-117, 123, 127-129, 140, 147-148, 154-156, 172-175, 186, 188, 189->191, 203, 222-226, 229-231, 235-238, 275-277, 297-311, 358, 382, 392-402, 408-413, 419, 427-441, 445-448, 454-456, 478->477, 501-578, 598-613, 624, 628, 646-651, 659-664, 668-670, 677-684 |
| src/plinder/data/structure/contacts.py                        |       88 |        8 |       30 |        5 |     88.9831% |93->95, 112-113, 191->195, 207, 220, 226-232, 258-261 |
| src/plinder/data/structure/models.py                          |       23 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/utils/annotations/\_\_init\_\_.py            |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/utils/annotations/aggregate\_annotations.py  |      505 |       73 |      339 |       63 |     79.1469% |85, 132->131, 142->141, 152->151, 162->161, 166->165, 170->169, 174->173, 175, 178->177, 182->181, 193->192, 196->exit, 199, 202->201, 206->205, 210->209, 219->218, 220-231, 234->233, 235-242, 256, 265, 273->272, 274, 277->276, 278, 320->319, 328->327, 336->335, 343->342, 347->346, 349, 357->356, 359, 384->390, 385->384, 387->385, 426-427, 466, 475, 481, 489-493, 502-503, 521, 532, 544, 610, 690->692, 703->702, 715->717, 725->724, 789->794, 792-793, 818->820, 822->826, 841->843, 867->879, 914->903, 934->933, 945, 978-999, 1002, 1007, 1010, 1149-1150, 1212-1215, 1218-1219, 1225->exit, 1230, 1277, 1300, 1304, 1316-1327 |
| src/plinder/data/utils/annotations/extras.py                  |      306 |       90 |       86 |       11 |     68.6224% |43->47, 66-69, 106-108, 123, 186->exit, 227, 259->256, 289, 307-312, 324, 344-348, 352-353, 357-365, 369-379, 406-415, 493->499, 521->526, 543, 556-557, 563-564, 572-573, 640->646, 655-698, 702-704, 710-725, 731-735 |
| src/plinder/data/utils/annotations/get\_ligand\_validation.py |      159 |       13 |       32 |       11 |     87.4346% |91->90, 92, 147->146, 154-158, 170, 197-201, 225->224, 234, 330->329, 338-339, 342, 361->363, 362->361, 363->362, 366 |
| src/plinder/data/utils/annotations/get\_similarity\_scores.py |      519 |      408 |      268 |       10 |     17.9161% |81-129, 133-136, 150, 207->220, 221, 227->230, 232, 233->230, 242-279, 312, 315->314, 322-323, 365-375, 383-391, 396-403, 412-497, 504-528, 533-596, 599-640, 643-650, 664-744, 755-809, 814-819, 824-900, 903-906, 909-918, 925-986, 991-1039 |
| src/plinder/data/utils/annotations/interaction\_utils.py      |      166 |       18 |       92 |        5 |     85.6589% |57-63, 70-75, 150-151, 199-202, 254, 363, 392->400 |
| src/plinder/data/utils/annotations/interface\_gap.py          |       87 |        8 |       28 |        6 |     87.8261% |53, 117, 126, 128, 164-165, 183-184 |
| src/plinder/data/utils/annotations/ligand\_utils.py           |      567 |       62 |      272 |       53 |     84.8629% |110->109, 133->exit, 165, 174, 189, 235-237, 260, 275-280, 302->291, 310, 332, 334, 336, 339-342, 346->345, 433->432, 439->441, 448->447, 458->457, 500-502, 546, 551, 577-578, 590-591, 602-603, 614-615, 627-628, 639-640, 713, 840->845, 847-850, 882-885, 918->917, 985, 1105-1106, 1108, 1128->1127, 1131->1133, 1136->1135, 1143->1142, 1147->1146, 1153->1152, 1159->1158, 1166->1165, 1167, 1173->1172, 1176->1178, 1182, 1204->1203, 1208->1207, 1215->1214, 1217, 1221->1220, 1225->1224, 1229->1228, 1233->1232, 1234-1241, 1244->1243, 1247-1248, 1252->1251, 1259-1260, 1293, 1328, 1362->1364, 1457->1461 |
| src/plinder/data/utils/annotations/mmpdb\_utils.py            |      124 |       66 |       40 |        0 |     45.1220% |35, 47-165, 173-185, 405-414 |
| src/plinder/data/utils/annotations/protein\_utils.py          |      137 |        8 |       70 |       12 |     90.3382% |47, 75->81, 79-80, 107->105, 159->161, 222, 238, 268->267, 269, 288, 336->335, 369->371, 382->381, 386->385, 387, 395->397, 415->414 |
| src/plinder/data/utils/annotations/rdkit\_utils.py            |      196 |       49 |       54 |        9 |     73.6000% |96-104, 146, 148, 153-158, 169->172, 174->161, 202-208, 212-219, 247-248, 289-296, 318-319, 369-370, 383-392, 397, 399-403, 414, 417->426, 422-423 |
| src/plinder/data/utils/annotations/save\_utils.py             |       80 |        4 |       38 |        5 |     90.6780% |39->28, 88-90, 115->117, 118->exit, 148 |
| src/plinder/data/utils/cluster.py                             |       43 |       43 |       12 |        0 |      0.0000% |     3-101 |
| src/plinder/data/utils/diffdock\_utils.py                     |      177 |      147 |       60 |        0 |     12.6582% |20, 26, 32, 38-43, 58-79, 86-93, 96-98, 102-121, 127-169, 173-188, 372-375, 379-446, 450-464, 476-520 |
| src/plinder/data/utils/tanimoto.py                            |      105 |       11 |       26 |        8 |     85.4962% |19, 31-32, 56-57, 62, 75, 137, 191, 193, 195 |
| src/plinder/eval/\_\_init\_\_.py                              |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/eval/docking/\_\_init\_\_.py                      |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/eval/docking/make\_plots.py                       |      113 |      113 |       24 |        0 |      0.0000% |     3-265 |
| src/plinder/eval/docking/stratify\_test\_set.py               |      117 |      117 |       38 |        0 |      0.0000% |     3-372 |
| src/plinder/eval/docking/utils.py                             |      184 |        5 |       68 |       15 |     92.0635% |46->45, 106->105, 117, 123, 147->149, 255->244, 264->287, 287->290, 298, 303->301, 307, 317, 336->338, 343->345, 345->360 |
| src/plinder/eval/docking/write\_scores.py                     |       84 |       84 |       38 |        0 |      0.0000% |     3-224 |
| src/plinder/methods/\_\_init\_\_.py                           |        0 |        0 |        0 |        0 |    100.0000% |           |
|                                                     **TOTAL** | **7973** | **2899** | **3095** |  **495** | **59.9837%** |           |


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