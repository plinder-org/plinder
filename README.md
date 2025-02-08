# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/plinder-org/plinder/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                                           |    Stmts |     Miss |   Branch |   BrPart |        Cover |   Missing |
|--------------------------------------------------------------- | -------: | -------: | -------: | -------: | -----------: | --------: |
| src//plinder/\_\_init\_\_.py                                   |        4 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/\_version.py                                      |       15 |        4 |        0 |        0 |     73.3333% |     16-19 |
| src//plinder/core/\_\_init\_\_.py                              |        5 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/index/\_\_init\_\_.py                        |        2 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/index/system.py                              |      198 |       55 |       50 |        6 |     66.5323% |13, 40, 75-76, 112, 211-212, 257-260, 272-281, 299-323, 337-338, 352-353, 391-405, 425, 430->435, 431->430, 436 |
| src//plinder/core/index/utils.py                               |      133 |       18 |       58 |       11 |     82.7225% |75, 162, 170-175, 177-178, 180-181, 199, 252, 260, 275, 291, 293 |
| src//plinder/core/loader/\_\_init\_\_.py                       |        2 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/loader/dataset.py                            |       35 |        3 |        4 |        2 |     87.1795% |59, 65, 70->73, 94 |
| src//plinder/core/loader/featurizer.py                         |       29 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/loader/transforms.py                         |       13 |        2 |        0 |        0 |     84.6154% |     6, 12 |
| src//plinder/core/loader/utils.py                              |       56 |       29 |       12 |        3 |     44.1176% |41->44, 48, 98->102, 116-128, 178-210 |
| src//plinder/core/scores/\_\_init\_\_.py                       |        8 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/scores/clusters.py                           |       19 |        2 |        2 |        1 |     85.7143% |     51-52 |
| src//plinder/core/scores/index.py                              |       29 |        3 |       10 |        3 |     84.6154% |42, 46, 60 |
| src//plinder/core/scores/ligand.py                             |       46 |        3 |        8 |        2 |     90.7407% | 49-52, 66 |
| src//plinder/core/scores/links.py                              |       32 |        7 |       10 |        4 |     73.8095% |44-45, 47-48, 51, 65-66 |
| src//plinder/core/scores/protein.py                            |       77 |        6 |       30 |        8 |     86.9159% |47->53, 61-62, 160, 164->174, 170, 173, 181->183, 190 |
| src//plinder/core/scores/query.py                              |      104 |        9 |       56 |        8 |     88.1250% |49, 97, 104, 173->181, 176-178, 216, 218, 264 |
| src//plinder/core/split/\_\_init\_\_.py                        |        2 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/split/plot.py                                |      302 |      182 |       76 |       10 |     35.4497% |19-20, 39, 156, 158, 159->161, 171, 176, 189-205, 295-313, 329-330, 379-397, 400-428, 431-473, 476-603, 606-614, 625-626, 629-660, 663-710, 713-755, 762, 764-766, 769-773, 777-820, 832 |
| src//plinder/core/split/utils.py                               |       19 |        8 |        2 |        0 |     52.3810% |     35-42 |
| src//plinder/core/structure/\_\_init\_\_.py                    |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/structure/atoms.py                           |      178 |       43 |       56 |        9 |     69.2308% |85, 99-101, 115, 139-140, 159-162, 166, 218-248, 283->282, 292-296, 300-302, 306, 335, 389-396 |
| src//plinder/core/structure/contacts.py                        |       88 |        8 |       30 |        5 |     88.9831% |93->95, 112-113, 191->195, 207, 220, 226-232, 258-261 |
| src//plinder/core/structure/diffdock\_utils.py                 |      177 |      137 |       58 |        1 |     18.2979% |26, 32, 38-43, 58-79, 86-93, 96-98, 102-121, 127-169, 173-188, 374-375, 388, 450-464, 476-520 |
| src//plinder/core/structure/models.py                          |       12 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/structure/structure.py                       |      328 |      129 |       80 |       16 |     55.6373% |37, 63-73, 82, 86-89, 165, 175, 181-184, 201-209, 213, 215, 217, 242-243, 246-268, 285-287, 298, 300->303, 315-317, 323-326, 336-399, 402-403, 415, 417, 470-473, 529, 548-553, 558-563, 598-599, 623->628, 633-635, 640-642, 647-649, 654-660, 665-668, 673-679, 684-688, 693-697, 711-715, 722->724, 724->726, 730 |
| src//plinder/core/structure/superimpose.py                     |       57 |        7 |       16 |        5 |     80.8219% |92, 104-107, 115, 138, 162->166 |
| src//plinder/core/structure/surgery.py                         |       65 |       51 |       18 |        0 |     16.8675% |19-30, 38-48, 54-86, 97, 118-129 |
| src//plinder/core/structure/vendored.py                        |      367 |      213 |      116 |       17 |     39.1304% |30-32, 36-38, 44-53, 70-72, 77-87, 96, 103-104, 124-127, 137->139, 140, 141->143, 161-162, 167-173, 179-184, 194-284, 342, 352, 377->379, 379->381, 384, 424->427, 440->442, 442->444, 444->439, 483-503, 515, 517, 531-551, 558-578, 613-638, 671-688, 694-707, 763-798, 817-826, 842-844, 862-865, 885-903 |
| src//plinder/core/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/utils/config.py                              |      125 |        4 |       40 |        5 |     94.5455% |130->129, 134, 155, 262, 284 |
| src//plinder/core/utils/constants.py                           |       48 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/utils/cpl.py                                 |      109 |       37 |       34 |        7 |     59.4406% |29-37, 46, 51, 58, 69, 74-85, 102-105, 114, 152->154, 159-168, 178 |
| src//plinder/core/utils/dataclass.py                           |       33 |        3 |       18 |        4 |     86.2745% |37, 47, 54->57, 65 |
| src//plinder/core/utils/dec.py                                 |       19 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/utils/gcs.py                                 |       93 |       25 |       24 |        7 |     67.5214% |28->50, 33-35, 47, 54-58, 117-120, 132-133, 135, 145, 157-164 |
| src//plinder/core/utils/io.py                                  |       47 |        7 |        6 |        2 |     83.0189% |54->62, 57-58, 67-81 |
| src//plinder/core/utils/load\_systems.py                       |       23 |       23 |       12 |        0 |      0.0000% |      4-41 |
| src//plinder/core/utils/log.py                                 |       37 |       12 |       14 |        2 |     64.7059% |12-13, 63-66, 83-88 |
| src//plinder/core/utils/schemas.py                             |        9 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/core/utils/unpack.py                              |       88 |        8 |       44 |        7 |     88.6364% |33, 39, 45, 114, 127->132, 134->133, 138-140, 150 |
| src//plinder/data/\_\_init\_\_.py                              |        6 |        2 |        0 |        0 |     66.6667% |       8-9 |
| src//plinder/data/\_version.py                                 |       15 |       15 |        0 |        0 |      0.0000% |     18-40 |
| src//plinder/data/clusters.py                                  |      118 |       37 |       26 |        6 |     65.9722% |195-232, 257-258, 267->exit, 274, 305-306, 320-321, 332-369, 380 |
| src//plinder/data/column\_descriptions/\_\_init\_\_.py         |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/common/\_\_init\_\_.py                       |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/common/\_version.py                          |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/common/constants.py                          |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/common/log.py                                |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/databases.py                                 |       59 |        3 |       20 |        2 |     93.6709% |13, 137->148, 153-154 |
| src//plinder/data/docs.py                                      |       57 |        0 |       18 |        1 |     98.6667% |    85->87 |
| src//plinder/data/final\_structure\_qc.py                      |      152 |       34 |       28 |        6 |     76.6667% |22, 45-47, 69-72, 114-123, 144-145, 147, 151-152, 154, 252-253, 280-283, 367-372 |
| src//plinder/data/get\_system\_annotations.py                  |      106 |       46 |       34 |        5 |     49.2857% |59, 62-63, 73-112, 141->145, 154->exit, 158-175, 187-188 |
| src//plinder/data/leakage.py                                   |       93 |       83 |       14 |        0 |      9.3458% |14-15, 27-51, 59-90, 100-163, 177-229 |
| src//plinder/data/pipeline/\_\_init\_\_.py                     |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/pipeline/config.py                           |      117 |        3 |       22 |        3 |     95.6835% |160, 183, 187 |
| src//plinder/data/pipeline/io.py                               |      272 |       60 |       78 |       11 |     74.5714% |119-145, 269->277, 363->361, 433, 439, 444, 448, 469, 474-476, 500-505, 515-524, 547-576 |
| src//plinder/data/pipeline/mpqueue.py                          |       65 |       15 |       12 |        5 |     74.0260% |27, 32-35, 51-56, 63->65, 80-85, 88->exit, 135->134 |
| src//plinder/data/pipeline/pipeline.py                         |      181 |       39 |       14 |        1 |     78.4615% |52-56, 193, 246-249, 258-259, 263, 267-279, 285-288, 296-300, 304, 308-312, 316, 323-329, 333, 340-344, 348, 355-359, 368-373, 377-388, 392-393, 439-440 |
| src//plinder/data/pipeline/tasks.py                            |      379 |      151 |      110 |       12 |     55.6237% |161, 174->177, 178-179, 218->227, 287-288, 297-298, 327-328, 334, 341, 372, 403->405, 409-413, 414->416, 424-427, 428->395, 466, 723-731, 737-738, 753-772, 793-812, 821-822, 844-864, 874-914, 922-923, 931-947, 955-977, 985-986, 999-1000, 1013, 1022-1028, 1044-1070, 1079-1090, 1101-1125 |
| src//plinder/data/pipeline/transform.py                        |      140 |       81 |       64 |        4 |     41.6667% |26, 29-30, 34, 48, 126-170, 188-211, 235-320, 357 |
| src//plinder/data/pipeline/utils.py                            |      375 |      134 |      128 |       13 |     62.0278% |28, 50-53, 103-104, 142-143, 155-156, 165, 167->170, 171-172, 174-175, 206-207, 229->227, 293, 324-331, 479-500, 509, 568->565, 590, 594->596, 605-610, 621-651, 662-686, 705-736, 751-752, 759-767, 780-797, 811-820 |
| src//plinder/data/save\_linked\_structures.py                  |      163 |       99 |       44 |        4 |     33.8164% |23-32, 37-44, 70-108, 136-137, 148, 180->182, 187-188, 223-224, 236-246, 258-290, 301-327, 340-347, 383-385 |
| src//plinder/data/smallmolecules.py                            |      110 |       36 |       14 |        6 |     64.5161% |19-20, 26-27, 33-34, 40-41, 47-48, 54-55, 61-62, 71-72, 83-89, 101, 105-106, 116-118, 125->124, 127-129, 139-140, 152->154, 163-166 |
| src//plinder/data/splits.py                                    |      298 |      233 |       58 |        0 |     18.2584% |121-126, 144-151, 155-162, 166-169, 191-199, 210-222, 243-334, 355-394, 418-468, 501-519, 561-619, 656-901, 911-970 |
| src//plinder/data/structure/\_\_init\_\_.py                    |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/structure/atoms.py                           |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/structure/contacts.py                        |        3 |        3 |        0 |        0 |      0.0000% |       3-7 |
| src//plinder/data/utils/\_\_init\_\_.py                        |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/utils/annotations/\_\_init\_\_.py            |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/data/utils/annotations/aggregate\_annotations.py  |      593 |       79 |      244 |       35 |     82.7957% |237, 270-281, 288-295, 307, 315, 492, 505, 539, 552, 565, 593->591, 632-633, 664, 670, 678-682, 691-692, 710, 721, 733, 833->835, 934->939, 937-938, 971->973, 975->979, 994->996, 1020->1032, 1067->1056, 1091, 1124-1138, 1141, 1146, 1149, 1184, 1276-1277, 1339-1342, 1345-1346, 1352->exit, 1395, 1418, 1422, 1434-1445, 1454-1475 |
| src//plinder/data/utils/annotations/get\_ligand\_validation.py |      137 |       12 |       18 |        4 |     89.6774% |66-70, 82, 123-127, 196, 307-308, 311, 338 |
| src//plinder/data/utils/annotations/get\_similarity\_scores.py |      519 |      408 |      234 |        9 |     17.5299% |81-129, 133-136, 150, 207->220, 221, 227->230, 232, 233->230, 242-279, 312, 322-323, 365-375, 383-391, 396-403, 412-497, 504-528, 533-596, 599-640, 643-650, 664-744, 755-809, 814-819, 824-900, 903-906, 909-918, 925-986, 991-1039 |
| src//plinder/data/utils/annotations/interaction\_utils.py      |      148 |        6 |       80 |        5 |     93.4211% |109-110, 158-161, 213, 322, 351->359 |
| src//plinder/data/utils/annotations/interface\_gap.py          |       87 |        8 |       28 |        6 |     87.8261% |53, 119, 128, 130, 166-167, 185-186 |
| src//plinder/data/utils/annotations/ligand\_utils.py           |      588 |       59 |      204 |       27 |     88.1313% |175, 184, 199, 245-247, 285-290, 312->301, 320, 342, 344, 346, 349-352, 516-518, 562, 567, 593-594, 606-607, 618-619, 630-631, 643-644, 655-656, 886->891, 893-896, 928-931, 1017, 1141-1142, 1144, 1170->1172, 1237->1239, 1243, 1287, 1316-1323, 1332-1333, 1344-1345, 1378, 1435, 1472->1474, 1527->1531 |
| src//plinder/data/utils/annotations/mmpdb\_utils.py            |      124 |       66 |       12 |        0 |     42.6471% |35, 47-165, 173-185, 405-414 |
| src//plinder/data/utils/annotations/protein\_utils.py          |      138 |        4 |       52 |        8 |     93.6842% |61-62, 88->86, 95->102, 97->95, 99->95, 153->155, 266, 337->339, 361, 388->387 |
| src//plinder/data/utils/annotations/rdkit\_utils.py            |      265 |       78 |       78 |       13 |     68.2216% |53-61, 103, 105, 110-115, 126->129, 131->118, 159-165, 169-176, 204-205, 288-289, 313->346, 322->337, 338-343, 349, 362-369, 391-392, 442-443, 457-466, 471, 473-477, 488, 491->500, 496-497, 505-532 |
| src//plinder/data/utils/annotations/save\_utils.py             |       80 |        4 |       30 |        2 |     92.7273% |88-90, 148 |
| src//plinder/data/utils/annotations/utils.py                   |       42 |        2 |       22 |        3 |     92.1875% |31->33, 51, 53 |
| src//plinder/data/utils/tanimoto.py                            |      105 |       11 |       26 |        8 |     85.4962% |19, 31-32, 56-57, 62, 75, 137, 192, 194, 196 |
| src//plinder/eval/\_\_init\_\_.py                              |        5 |        2 |        0 |        0 |     60.0000% |       7-8 |
| src//plinder/eval/docking/\_\_init\_\_.py                      |        0 |        0 |        0 |        0 |    100.0000% |           |
| src//plinder/eval/docking/make\_plots.py                       |      126 |        5 |       26 |        5 |     93.4211% |76->79, 85-86, 120, 265, 276 |
| src//plinder/eval/docking/stratify\_test\_set.py               |      111 |       11 |       28 |        8 |     84.8921% |68->75, 141, 205->204, 249-252, 254-257, 290-293, 345, 358 |
| src//plinder/eval/docking/utils.py                             |      226 |       19 |       70 |       14 |     88.1757% |81, 226-242, 341->330, 348->394, 384-390, 394->399, 407, 412->410, 416, 427, 442, 446-448, 453, 482->485, 494->500, 501->486 |
| src//plinder/eval/docking/write\_scores.py                     |      107 |       53 |       44 |        4 |     42.3841% |53-76, 116-167, 209->214, 238-239, 243->242, 303, 321 |
| src//plinder/methods/\_\_init\_\_.py                           |        0 |        0 |        0 |        0 |    100.0000% |           |
|                                                      **TOTAL** | **9152** | **2948** | **2764** |  **375** | **64.4092%** |           |


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