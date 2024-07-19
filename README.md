# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/plinder-org/plinder/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                                          |    Stmts |     Miss |   Branch |   BrPart |        Cover |   Missing |
|-------------------------------------------------------------- | -------: | -------: | -------: | -------: | -----------: | --------: |
| src/plinder/core/\_\_init\_\_.py                              |        4 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/scores/\_\_init\_\_.py                       |        4 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/scores/clusters.py                           |       18 |        7 |        4 |        1 |     54.5455% |18->17, 39-52 |
| src/plinder/core/scores/ligand.py                             |       18 |        0 |        4 |        1 |     95.4545% |    18->17 |
| src/plinder/core/scores/protein.py                            |       25 |        0 |       12 |        2 |     94.5946% |18->17, 45->50 |
| src/plinder/core/scores/query.py                              |       85 |       34 |       38 |        9 |     53.6585% |36-37, 65, 86-88, 129, 133, 136, 142, 148, 155->157, 158, 191-219 |
| src/plinder/core/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/core/utils/config.py                              |      110 |        6 |       40 |        3 |     92.6667% |122, 142, 146-147, 267-268 |
| src/plinder/core/utils/cpl.py                                 |       63 |       17 |       26 |       12 |     65.1685% |21->20, 26-30, 38->37, 39, 43->42, 44, 51, 55->57, 56->55, 57->56, 62->61, 72-79, 102, 104, 106 |
| src/plinder/core/utils/dec.py                                 |       19 |        0 |        2 |        1 |     95.2381% |    18->17 |
| src/plinder/core/utils/gcs.py                                 |       93 |       25 |       40 |       14 |     64.6617% |26->25, 28->50, 33-35, 47, 54-58, 66->65, 79->78, 93->92, 109->108, 117-120, 124->123, 132-133, 135, 142->exit, 145, 157-164 |
| src/plinder/core/utils/log.py                                 |       37 |       12 |       14 |        2 |     64.7059% |12-13, 63-66, 83-88 |
| src/plinder/core/utils/schemas.py                             |        9 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/\_\_init\_\_.py                              |        4 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/\_version.py                                 |       15 |        4 |        0 |        0 |     73.3333% |     35-39 |
| src/plinder/data/clusters.py                                  |       80 |       63 |       16 |        0 |     17.7083% |38-44, 68-73, 110-135, 167-186, 201-232, 242-301 |
| src/plinder/data/common/\_\_init\_\_.py                       |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/common/constants.py                          |       30 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/common/log.py                                |        3 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/databases.py                                 |       59 |       45 |       26 |        1 |     17.6471% |13, 19-20, 49-60, 93-149, 153-154, 181-187, 209-222 |
| src/plinder/data/final\_structure\_qc.py                      |      144 |       30 |       28 |        7 |     77.3256% |22, 45-47, 69-72, 114-123, 147, 245-246, 273-276, 360-365, 457, 480->482 |
| src/plinder/data/get\_system\_annotations.py                  |      105 |       71 |       44 |        4 |     25.5034% |56, 66, 69-70, 80-120, 124-181, 185-202 |
| src/plinder/data/leakage.py                                   |       93 |       83 |       16 |        0 |      9.1743% |14-15, 27-51, 59-90, 100-163, 177-229 |
| src/plinder/data/pipeline/\_\_init\_\_.py                     |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/pipeline/config.py                           |      120 |        4 |       44 |        7 |     93.2927% |79->83, 83->exit, 127, 146->148, 151, 155, 195, 257->252 |
| src/plinder/data/pipeline/io.py                               |      311 |       76 |      136 |       34 |     67.7852% |40->39, 60->59, 90->95, 93->95, 99->98, 147-194, 196->198, 202->201, 234->237, 235->234, 246->245, 277->281, 285->284, 318->326, 330->329, 359->361, 365->364, 411->409, 461->460, 481, 486-487, 491, 505->504, 512, 517-519, 543-548, 558-567, 571->570, 578->586, 581-582, 584->586, 590->589, 609-638, 642->641, 670->672 |
| src/plinder/data/pipeline/pipeline.py                         |      152 |       60 |       78 |       33 |     57.8261% |52-61, 64->63, 65-70, 73->72, 74, 80->79, 81, 88->87, 89, 95->94, 96-106, 109->108, 110-121, 124->123, 125-128, 131->130, 132-137, 140->139, 141, 147->146, 148-158, 161->160, 162, 168->167, 169, 175->174, 176-181, 184->183, 185, 194->193, 195, 198->197, 199-204, 207->206, 208, 214->213, 215, 221->220, 222-228, 231->230, 232, 239->238, 240-246, 249->248, 250, 257->256, 258, 261->260, 262-268, 271->270, 274, 281->280, 282-286, 289->288, 290, 293->292, 294-298, 301->300, 302, 308->307, 309-315, 318->317, 319, 325->324, 326, 376-377 |
| src/plinder/data/pipeline/tasks.py                            |      311 |      250 |       84 |        0 |     16.4557% |139-153, 165-175, 213-232, 284-365, 391-397, 405-454, 485-491, 504-516, 536-541, 556-580, 591-597, 611-615, 630-635, 648-651, 688-708, 740-760, 772-788, 797-805, 826-833, 850-853, 862-863, 885-905, 915-955, 963-964, 972-988, 996-1018, 1026-1027, 1040-1041, 1055-1059 |
| src/plinder/data/pipeline/transform.py                        |      140 |       81 |       85 |        5 |     38.6667% |26, 29-30, 34, 48, 99->101, 125-168, 186-209, 233-318, 355 |
| src/plinder/data/pipeline/utils.py                            |      232 |       95 |      110 |       11 |     56.1404% |22-23, 36->35, 45-48, 83->82, 109-110, 116->115, 148-149, 150->146, 154->150, 156->154, 161-162, 170-182, 192-218, 233-244, 296-298, 326-333, 380->379, 395-411, 420-423, 432-459 |
| src/plinder/data/smallmolecules.py                            |      101 |       34 |       12 |        5 |     63.7168% |19-20, 26-27, 33-34, 40-41, 47-48, 54-55, 61-62, 71-72, 83-89, 101, 105-106, 116-118, 125->124, 127-129, 139-140 |
| src/plinder/data/splits.py                                    |      293 |      215 |       86 |        2 |     22.6913% |56->exit, 57->exit, 108-113, 131-138, 142-150, 154-157, 161-195, 217-225, 236-250, 269-362, 383-423, 448-504, 538-556, 604-653, 690-761, 765-820 |
| src/plinder/data/structure/\_\_init\_\_.py                    |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/structure/atoms.py                           |      301 |      148 |       78 |       20 |     47.2296% |73-75, 85-96, 103, 109->113, 115-117, 123, 127-129, 140, 147-148, 154-156, 172-175, 186, 188, 189->191, 203, 222-226, 229-231, 235-238, 275-277, 297-311, 358, 382, 392-402, 408-413, 419, 427-441, 445-448, 454-456, 478->477, 501-578, 598-613, 624, 628, 646-651, 659-664, 668-670, 677-684 |
| src/plinder/data/structure/contacts.py                        |       88 |        8 |       30 |        5 |     88.9831% |93->95, 112-113, 191->195, 207, 220, 226-232, 258-261 |
| src/plinder/data/structure/models.py                          |       23 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/utils/annotations/\_\_init\_\_.py            |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/data/utils/annotations/aggregate\_annotations.py  |      463 |       77 |      294 |       49 |     77.0145% |79, 126->125, 136->135, 146->145, 150->149, 154->153, 158->157, 169->168, 178->177, 182->181, 191->190, 192-201, 204->203, 205-212, 226, 235, 243->242, 244, 247->246, 248, 283->282, 291->290, 306->312, 307->306, 309->307, 351->353, 364-365, 400, 409, 415, 423-427, 436-439, 457, 468, 480, 545, 625->627, 638->637, 650->652, 660->659, 722->727, 725-726, 750->752, 754->758, 783->785, 809->821, 856->845, 872->866, 876->875, 887, 908-952, 1091-1092, 1154-1155, 1158-1161, 1167->exit, 1172, 1202, 1225, 1229, 1241-1252 |
| src/plinder/data/utils/annotations/extras.py                  |      397 |      124 |      118 |       15 |     66.0194% |44->48, 67-70, 107-109, 124, 187->exit, 228, 260->257, 290, 308-313, 325, 345-349, 353-354, 358-366, 370-380, 407-416, 494->500, 522->527, 544, 557-558, 564-565, 573-574, 641->647, 654-709, 726->723, 736->735, 795-798, 833, 870-913, 917-919, 925-940, 946-950 |
| src/plinder/data/utils/annotations/get\_ligand\_validation.py |      153 |       14 |       32 |       12 |     85.9459% |75->74, 76, 120->119, 127-131, 143, 170-174, 198->197, 207, 293->292, 298, 301-302, 305, 324->326, 325->324, 326->325, 329 |
| src/plinder/data/utils/annotations/interaction\_utils.py      |      180 |       24 |       86 |        4 |     84.2105% |57-63, 70-75, 125-131, 222-226, 277, 395, 424->432 |
| src/plinder/data/utils/annotations/interface\_gap.py          |       87 |        8 |       28 |        6 |     87.8261% |53, 117, 126, 128, 164-165, 183-184 |
| src/plinder/data/utils/annotations/ligand\_utils.py           |      515 |      107 |      238 |       47 |     76.3612% |98->97, 106-125, 140-146, 153, 162, 177, 223-225, 248, 263-268, 290->279, 298, 323, 325, 327, 330-333, 337->336, 346-420, 424->423, 429-435, 439->438, 440-445, 449->448, 491-493, 537, 542, 568-569, 581-582, 593-594, 605-606, 618-619, 630-631, 693, 799->805, 807-810, 842-845, 878->877, 940-944, 946, 948, 950, 1058-1059, 1061, 1079->1078, 1087->1086, 1091->1090, 1097->1096, 1103->1102, 1110->1109, 1113->1115, 1119, 1125->1124, 1129->1128, 1133->1132, 1137->1136, 1138-1145, 1148->1147, 1151-1154, 1158->1157, 1165-1168, 1198, 1233, 1267->1269, 1356->1360 |
| src/plinder/data/utils/annotations/mmpdb\_utils.py            |      124 |       66 |       40 |        0 |     45.1220% |35, 47-165, 173-185, 405-414 |
| src/plinder/data/utils/annotations/protein\_utils.py          |      137 |        8 |       70 |       12 |     90.3382% |47, 75->81, 79-80, 107->105, 159->161, 222, 238, 268->267, 269, 288, 336->335, 369->371, 382->381, 386->385, 387, 395->397, 415->414 |
| src/plinder/data/utils/annotations/rdkit\_utils.py            |      182 |       35 |       48 |        6 |     78.6957% |92-100, 142, 144, 149-154, 165->168, 170->157, 198-204, 208-215, 243-244, 313-314, 352-353, 390->399, 395-396 |
| src/plinder/data/utils/annotations/save\_utils.py             |       80 |        1 |       38 |        5 |     94.9153% |35->27, 84->87, 111->113, 114->exit, 144 |
| src/plinder/data/utils/diffdock\_utils.py                     |      177 |      147 |       60 |        0 |     12.6582% |20, 26, 32, 38-43, 58-79, 86-93, 96-98, 102-121, 127-169, 173-188, 372-375, 379-446, 450-464, 476-520 |
| src/plinder/data/utils/tanimoto.py                            |      105 |       87 |       26 |        1 |     14.5038% |19, 27-32, 53-97, 103-171, 183-226 |
| src/plinder/eval/docking/\_\_init\_\_.py                      |        0 |        0 |        0 |        0 |    100.0000% |           |
| src/plinder/eval/docking/utils.py                             |      171 |        3 |       64 |       13 |     93.1915% |40->39, 100->99, 130->132, 238->227, 247->264, 264->267, 275, 280->278, 284, 294, 313->315, 320->322, 322->337 |
|                                                     **TOTAL** | **5861** | **2069** | **2195** |  **349** | **61.4449%** |           |


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