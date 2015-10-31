{
	"patcher" : 	{
		"fileversion" : 1,
		"rect" : [ 335.0, 44.0, 758.0, 652.0 ],
		"bglocked" : 0,
		"defrect" : [ 335.0, 44.0, 758.0, 652.0 ],
		"openrect" : [ 0.0, 0.0, 0.0, 0.0 ],
		"openinpresentation" : 1,
		"default_fontsize" : 12.0,
		"default_fontface" : 0,
		"default_fontname" : "Arial",
		"gridonopen" : 0,
		"gridsize" : [ 5.0, 5.0 ],
		"gridsnaponopen" : 0,
		"toolbarvisible" : 1,
		"boxanimatetime" : 200,
		"imprint" : 0,
		"boxes" : [ 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-58",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 625.0, 465.0, 53.0, 20.0 ],
					"text" : "bitsafe~"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-24",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "float" ],
					"patching_rect" : [ 1280.0, 170.0, 34.0, 20.0 ],
					"text" : "/ 12."
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-16",
					"maxclass" : "flonum",
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "float", "bang" ],
					"patching_rect" : [ 1230.0, 170.0, 50.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 150.0, 300.0, 40.0, 20.0 ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-5",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 435.0, 180.0, 79.0, 20.0 ],
					"text" : "loadmess 40"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-80",
					"linecount" : 4,
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 10.0, 65.0, 150.0, 62.0 ],
					"presentation" : 1,
					"presentation_linecount" : 4,
					"presentation_rect" : [ 425.0, 75.0, 150.0, 62.0 ],
					"text" : "This is an LPC example that resynthesizes the input using a band-limited impulse train."
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-78",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 455.0, 50.0, 150.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 10.0, 120.0, 92.0, 20.0 ],
					"text" : "Soundfile Input:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial Bold",
					"fontsize" : 14.0,
					"id" : "obj-69",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 550.0, 15.0, 183.0, 23.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 10.0, 65.0, 183.0, 23.0 ],
					"text" : "Input Parameters:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-66",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 90.0, 165.0, 150.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 95.0, 200.0, 50.0, 20.0 ],
					"text" : "Preset:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial Bold",
					"fontsize" : 14.0,
					"id" : "obj-61",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 15.0, 320.0, 185.0, 23.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 15.0, 415.0, 185.0, 23.0 ],
					"text" : "Output:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial Bold",
					"fontsize" : 14.0,
					"id" : "obj-28",
					"linecount" : 2,
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 1065.0, 0.0, 183.0, 39.0 ],
					"presentation" : 1,
					"presentation_linecount" : 2,
					"presentation_rect" : [ 10.0, 150.0, 183.0, 39.0 ],
					"text" : "Analysis/Synthesis Parameters:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-23",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 1175.0, 775.0, 150.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 100.0, 360.0, 62.0, 20.0 ],
					"text" : "Listen:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-11",
					"maxclass" : "textbutton",
					"mode" : 1,
					"numinlets" : 1,
					"numoutlets" : 3,
					"outlettype" : [ "", "", "int" ],
					"patching_rect" : [ 550.0, 50.0, 100.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 205.0, 120.0, 100.0, 20.0 ],
					"text" : "play",
					"texton" : "stop"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-7",
					"maxclass" : "textbutton",
					"numinlets" : 1,
					"numoutlets" : 3,
					"outlettype" : [ "", "", "int" ],
					"patching_rect" : [ 655.0, 50.0, 100.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 105.0, 120.0, 100.0, 20.0 ],
					"text" : "choose input file"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-3",
					"maxclass" : "slider",
					"min" : 1.0,
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 515.0, 180.0, 190.0, 18.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 190.0, 240.0, 192.0, 20.0 ],
					"size" : 100.0
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-77",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 465.0, 85.0, 150.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 65.0, 95.0, 44.0, 20.0 ],
					"text" : "Input:"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-75",
					"maxclass" : "preset",
					"numinlets" : 1,
					"numoutlets" : 4,
					"outlettype" : [ "preset", "int", "preset", "int" ],
					"patching_rect" : [ 915.0, 50.0, 100.0, 40.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 150.0, 190.0, 100.0, 40.0 ],
					"preset_data" : [ 						{
							"number" : 1,
							"data" : [ 5, "obj-72", "slider", "float", 60.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 1, 5, "obj-3", "slider", "float", 39.0 ]
						}
, 						{
							"number" : 2,
							"data" : [ 5, "obj-72", "slider", "float", 60.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 3, 5, "obj-3", "slider", "float", 39.0 ]
						}
, 						{
							"number" : 3,
							"data" : [ 5, "obj-72", "slider", "float", 60.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 2, 5, "obj-3", "slider", "float", 39.0 ]
						}
, 						{
							"number" : 4,
							"data" : [ 5, "obj-72", "slider", "float", 60.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 39.0 ]
						}
, 						{
							"number" : 5,
							"data" : [ 5, "obj-72", "slider", "float", 99.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 39.0 ]
						}
, 						{
							"number" : 6,
							"data" : [ 5, "obj-72", "slider", "float", 99.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 92.0 ]
						}
, 						{
							"number" : 7,
							"data" : [ 5, "obj-72", "slider", "float", 0.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 92.0 ]
						}
, 						{
							"number" : 8,
							"data" : [ 5, "obj-72", "slider", "float", 0.0, 5, "obj-64", "slider", "float", 240.0, 5, "obj-9", "slider", "float", 249.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 58.0 ]
						}
, 						{
							"number" : 9,
							"data" : [ 5, "obj-72", "slider", "float", 0.0, 5, "obj-64", "slider", "float", 121.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 78.0 ]
						}
, 						{
							"number" : 10,
							"data" : [ 5, "obj-72", "slider", "float", 49.0, 5, "obj-64", "slider", "float", 121.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 78.0 ]
						}
, 						{
							"number" : 11,
							"data" : [ 5, "obj-72", "slider", "float", 61.0, 5, "obj-64", "slider", "float", 0.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 57.0 ]
						}
, 						{
							"number" : 12,
							"data" : [ 5, "obj-72", "slider", "float", 49.0, 5, "obj-64", "slider", "float", 0.0, 5, "obj-9", "slider", "float", 387.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 57.0 ]
						}
, 						{
							"number" : 13,
							"data" : [ 5, "obj-72", "slider", "float", 13.0, 5, "obj-64", "slider", "float", 480.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 7.0 ]
						}
, 						{
							"number" : 14,
							"data" : [ 5, "obj-72", "slider", "float", 13.0, 5, "obj-64", "slider", "float", 0.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 7.0 ]
						}
, 						{
							"number" : 15,
							"data" : [ 5, "obj-72", "slider", "float", 76.0, 5, "obj-64", "slider", "float", 0.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 7.0 ]
						}
, 						{
							"number" : 16,
							"data" : [ 5, "obj-72", "slider", "float", 76.0, 5, "obj-64", "slider", "float", 362.0, 5, "obj-9", "slider", "float", 261.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 3.0 ]
						}
, 						{
							"number" : 17,
							"data" : [ 5, "obj-72", "slider", "float", 48.0, 5, "obj-64", "slider", "float", 480.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 57.0 ]
						}
, 						{
							"number" : 18,
							"data" : [ 5, "obj-72", "slider", "float", 48.0, 5, "obj-64", "slider", "float", 362.0, 5, "obj-9", "slider", "float", 5.0, 5, "obj-22", "umenu", "int", 4, 5, "obj-3", "slider", "float", 57.0 ]
						}
 ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-74",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 1190.0, 215.0, 79.0, 20.0 ],
					"text" : "loadmess 60"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-73",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "float" ],
					"patching_rect" : [ 1190.0, 285.0, 41.0, 20.0 ],
					"text" : "/ 100."
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-71",
					"linecount" : 2,
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 1085.0, 245.0, 102.0, 34.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 25.0, 330.0, 119.0, 20.0 ],
					"text" : "Unvoiced Threshold:"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-72",
					"maxclass" : "slider",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 1190.0, 240.0, 192.0, 38.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 190.0, 330.0, 192.0, 20.0 ],
					"size" : 100.0
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-70",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 1190.0, 100.0, 85.0, 20.0 ],
					"text" : "loadmess 240"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-68",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 1190.0, 195.0, 97.0, 20.0 ],
					"text" : "expr pow(2\\,$f1)"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-67",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "float" ],
					"patching_rect" : [ 1190.0, 170.0, 34.0, 20.0 ],
					"text" : "/ 10."
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-63",
					"linecount" : 2,
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 1085.0, 130.0, 102.0, 34.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 15.0, 300.0, 134.0, 20.0 ],
					"text" : "Pitch Shift (semitones):"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-64",
					"maxclass" : "slider",
					"min" : -240.0,
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 1190.0, 125.0, 192.0, 38.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 190.0, 300.0, 192.0, 20.0 ],
					"size" : 481.0
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-62",
					"maxclass" : "flonum",
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "float", "bang" ],
					"patching_rect" : [ 1125.0, 475.0, 50.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 150.0, 330.0, 40.0, 20.0 ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-60",
					"items" : [ "none", ",", "sound", "file", ",", "adc" ],
					"maxclass" : "umenu",
					"numinlets" : 1,
					"numoutlets" : 3,
					"outlettype" : [ "int", "", "" ],
					"patching_rect" : [ 505.0, 85.0, 100.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 105.0, 95.0, 100.0, 20.0 ],
					"types" : [  ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-59",
					"linecount" : 2,
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 5.0, 535.0, 150.0, 34.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 335.0, 440.0, 234.0, 20.0 ],
					"text" : "LPC filter coefficient frequency response:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-57",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 10.0, 450.0, 62.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 14.0, 438.0, 62.0, 20.0 ],
					"text" : "spectrum:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-27",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 690.0, 585.0, 43.0, 20.0 ],
					"text" : "sig~ 1"
				}

			}
, 			{
				"box" : 				{
					"color" : [ 0.478431, 0.709804, 0.317647, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-26",
					"maxclass" : "newobj",
					"numinlets" : 4,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 585.0, 620.0, 125.0, 20.0 ],
					"text" : "mbc.allpole~ 40 10. 1"
				}

			}
, 			{
				"box" : 				{
					"color" : [ 0.478431, 0.709804, 0.317647, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-25",
					"maxclass" : "newobj",
					"numinlets" : 3,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 585.0, 510.0, 97.0, 20.0 ],
					"text" : "mbc.errfilt~ 40 1"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-22",
					"items" : [ "none", ",", "original", ",", "residual", "(error)", ",", "reconstructed", "(inverse", "filtered)", ",", "blit", "resynthesized" ],
					"maxclass" : "umenu",
					"numinlets" : 1,
					"numoutlets" : 3,
					"outlettype" : [ "int", "", "" ],
					"patching_rect" : [ 920.0, 620.0, 100.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 150.0, 360.0, 230.0, 20.0 ],
					"types" : [  ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-21",
					"maxclass" : "newobj",
					"numinlets" : 5,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 730.0, 685.0, 73.0, 20.0 ],
					"text" : "selector~ 4"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-20",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 435.0, 205.0, 76.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 75.0, 240.0, 76.0, 20.0 ],
					"text" : "LPC Order:"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-14",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 1085.0, 70.0, 102.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 50.0, 275.0, 102.0, 20.0 ],
					"text" : "Pitch Glide (ms):"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-9",
					"maxclass" : "slider",
					"mult" : 2.0,
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 1190.0, 60.0, 192.0, 38.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 190.0, 270.0, 192.0, 20.0 ],
					"size" : 500.0
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-17",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 295.0, 20.0, 130.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 610.0, 20.0, 130.0, 20.0 ],
					"text" : "markcartwright.com"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-18",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 295.0, 0.0, 137.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 610.0, 0.0, 137.0, 20.0 ],
					"text" : "by m. brozier cartwright"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial Bold",
					"fontsize" : 24.0,
					"id" : "obj-19",
					"maxclass" : "comment",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 13.0, 12.0, 201.0, 34.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 13.0, 12.0, 201.0, 34.0 ],
					"text" : "lpcToolkitDemo"
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 0.584314, 0.917647, 0.180392, 1.0 ],
					"id" : "obj-4",
					"maxclass" : "panel",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ -5.0, 0.0, 440.0, 60.0 ],
					"presentation" : 1,
					"presentation_rect" : [ -5.0, 0.0, 755.0, 60.0 ],
					"rounded" : 0
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-15",
					"maxclass" : "ezdac~",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 730.0, 795.0, 45.0, 45.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 665.0, 565.0, 65.0, 65.0 ]
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-13",
					"maxclass" : "gain~",
					"numinlets" : 2,
					"numoutlets" : 2,
					"outlettype" : [ "signal", "int" ],
					"patching_rect" : [ 730.0, 715.0, 46.0, 69.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 665.0, 460.0, 65.0, 104.0 ]
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-12",
					"maxclass" : "toggle",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "int" ],
					"patching_rect" : [ 630.0, 75.0, 20.0, 20.0 ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-10",
					"maxclass" : "message",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 655.0, 75.0, 37.0, 18.0 ],
					"text" : "open"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-8",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 2,
					"outlettype" : [ "signal", "bang" ],
					"patching_rect" : [ 640.0, 105.0, 49.0, 20.0 ],
					"save" : [ "#N", "sfplay~", "", 1, 120960, 0, "", ";" ],
					"text" : "sfplay~"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-6",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 690.0, 105.0, 47.0, 20.0 ],
					"text" : "adc~ 1"
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-2",
					"maxclass" : "meter~",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "float" ],
					"patching_rect" : [ 600.0, 155.0, 80.0, 13.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 205.0, 95.0, 100.0, 20.0 ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-1",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 1190.0, 35.0, 72.0, 20.0 ],
					"text" : "loadmess 5"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-29",
					"maxclass" : "newobj",
					"numinlets" : 3,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 585.0, 130.0, 126.0, 20.0 ],
					"text" : "selector~ 2"
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 0.866667, 0.866667, 0.866667, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"htextcolor" : [ 0.870588, 0.870588, 0.870588, 1.0 ],
					"id" : "obj-30",
					"maxclass" : "number",
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "int", "bang" ],
					"patching_rect" : [ 1220.0, 355.0, 37.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 150.0, 270.0, 40.0, 20.0 ],
					"triscale" : 0.9
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-31",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "float", "bang" ],
					"patching_rect" : [ 1200.0, 330.0, 39.0, 20.0 ],
					"text" : "t f b"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-32",
					"maxclass" : "newobj",
					"numinlets" : 3,
					"numoutlets" : 2,
					"outlettype" : [ "float", "float" ],
					"patching_rect" : [ 907.0, 301.0, 93.0, 20.0 ],
					"text" : "split 60. 800."
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 0.866667, 0.866667, 0.866667, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"htextcolor" : [ 0.870588, 0.870588, 0.870588, 1.0 ],
					"id" : "obj-33",
					"maxclass" : "flonum",
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "float", "bang" ],
					"patching_rect" : [ 1240.0, 385.0, 40.0, 20.0 ],
					"triscale" : 0.9
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-34",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1201.0, 409.0, 58.0, 20.0 ],
					"text" : "*~ 1."
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-35",
					"maxclass" : "button",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "bang" ],
					"patching_rect" : [ 1062.0, 556.0, 15.0, 15.0 ]
				}

			}
, 			{
				"box" : 				{
					"id" : "obj-36",
					"maxclass" : "button",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "bang" ],
					"patching_rect" : [ 1036.0, 556.0, 15.0, 15.0 ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-37",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "bang", "bang" ],
					"patching_rect" : [ 1036.0, 533.0, 45.0, 20.0 ],
					"text" : "edge~"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-38",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1149.0, 432.0, 52.0, 20.0 ],
					"text" : "*~ 0.1"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-39",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1029.0, 432.0, 36.0, 20.0 ],
					"text" : "sig~"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-40",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1184.0, 533.0, 59.0, 20.0 ],
					"text" : "==~ 0."
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-41",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1097.0, 505.0, 50.0, 20.0 ],
					"text" : "<~ 0.6"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-42",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1166.0, 556.0, 32.0, 20.0 ],
					"text" : "*~"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-43",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1097.0, 556.0, 37.0, 20.0 ],
					"text" : "*~"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-44",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1149.0, 410.0, 50.0, 20.0 ],
					"text" : "noise~"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-45",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 2,
					"outlettype" : [ "signal", "bang" ],
					"patching_rect" : [ 1201.0, 383.0, 39.0, 20.0 ],
					"text" : "line~"
				}

			}
, 			{
				"box" : 				{
					"color" : [ 0.478431, 0.709804, 0.317647, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-46",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "float", "float" ],
					"patching_rect" : [ 907.0, 266.0, 140.0, 20.0 ],
					"text" : "mbc.pitch~ 30 60 0.2"
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 1.0, 1.0, 1.0, 1.0 ],
					"bordercolor" : [ 0.0, 0.0, 0.0, 1.0 ],
					"fgcolor" : [ 0.878431, 0.878431, 0.878431, 1.0 ],
					"id" : "obj-47",
					"markercolor" : [ 0.721569, 0.721569, 0.721569, 1.0 ],
					"maxclass" : "spectroscope~",
					"numinlets" : 2,
					"numoutlets" : 0,
					"patching_rect" : [ 71.0, 432.0, 300.0, 100.0 ],
					"peakcolor" : [ 0.952941, 0.8, 0.8, 1.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 15.0, 460.0, 318.0, 171.0 ],
					"rounded" : 0,
					"sonomedcolor" : [ 0.0, 0.74902, 0.498039, 1.0 ],
					"sonomedhicolor" : [ 1.0, 0.74902, 0.0, 1.0 ],
					"sonomedlocolor" : [ 0.498039, 0.0, 0.498039, 1.0 ]
				}

			}
, 			{
				"box" : 				{
					"color" : [ 0.478431, 0.709804, 0.317647, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-48",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1166.0, 504.0, 95.0, 20.0 ],
					"text" : "mbc.blit~ 100"
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 0.819608, 0.819608, 0.819608, 1.0 ],
					"bwidthcolor" : [ 0.298039, 0.423529, 0.67451, 1.0 ],
					"fgcolor" : [ 0.666667, 0.666667, 0.666667, 1.0 ],
					"fontsize" : 12.0,
					"hbwidthcolor" : [ 0.823529, 0.290196, 0.211765, 1.0 ],
					"hcurvecolor" : [ 1.0, 0.086275, 0.086275, 1.0 ],
					"id" : "obj-49",
					"logfreq" : 0,
					"markercolor" : [ 0.509804, 0.509804, 0.509804, 1.0 ],
					"maxclass" : "filtergraph~",
					"nfilters" : 1,
					"numinlets" : 8,
					"numoutlets" : 7,
					"outlettype" : [ "list", "float", "float", "float", "float", "list", "int" ],
					"patching_rect" : [ 70.0, 533.0, 491.0, 265.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 335.0, 460.0, 330.0, 170.0 ],
					"setfilter" : [ 0, 0, 1, 0, 0, 40.0, 1.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-50",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 70.0, 387.0, 119.0, 20.0 ],
					"text" : "prepend highorder"
				}

			}
, 			{
				"box" : 				{
					"color" : [ 0.478431, 0.709804, 0.317647, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-51",
					"maxclass" : "newobj",
					"numinlets" : 3,
					"numoutlets" : 1,
					"outlettype" : [ "list" ],
					"patching_rect" : [ 70.0, 363.0, 158.0, 20.0 ],
					"text" : "mbc.coeffDisp~ 40 filter"
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 0.866667, 0.866667, 0.866667, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"htextcolor" : [ 0.870588, 0.870588, 0.870588, 1.0 ],
					"id" : "obj-52",
					"maxclass" : "number",
					"maximum" : 199,
					"minimum" : 1,
					"numinlets" : 1,
					"numoutlets" : 2,
					"outlettype" : [ "int", "bang" ],
					"patching_rect" : [ 512.0, 203.0, 37.0, 20.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 150.0, 240.0, 40.0, 20.0 ],
					"triscale" : 0.9
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-53",
					"maxclass" : "message",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "" ],
					"patching_rect" : [ 512.0, 224.0, 64.0, 18.0 ],
					"text" : "order $1"
				}

			}
, 			{
				"box" : 				{
					"color" : [ 0.478431, 0.709804, 0.317647, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-54",
					"maxclass" : "newobj",
					"numinlets" : 4,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 730.0, 620.0, 125.0, 20.0 ],
					"text" : "mbc.allpole~ 40 10. 1"
				}

			}
, 			{
				"box" : 				{
					"color" : [ 0.478431, 0.709804, 0.317647, 1.0 ],
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-55",
					"maxclass" : "newobj",
					"numinlets" : 1,
					"numoutlets" : 5,
					"outlettype" : [ "signal", "signal", "signal", "signal", "signal" ],
					"patching_rect" : [ 585.0, 260.0, 124.0, 20.0 ],
					"text" : "mbc.lpc~ 40 100 1"
				}

			}
, 			{
				"box" : 				{
					"fontname" : "Arial",
					"fontsize" : 12.0,
					"id" : "obj-56",
					"maxclass" : "newobj",
					"numinlets" : 2,
					"numoutlets" : 1,
					"outlettype" : [ "signal" ],
					"patching_rect" : [ 1097.0, 595.0, 52.0, 20.0 ],
					"text" : "*~ 1."
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 0.678431, 0.819608, 0.819608, 1.0 ],
					"id" : "obj-81",
					"maxclass" : "panel",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 870.0, 0.0, 557.0, 707.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 5.0, 145.0, 385.0, 260.0 ],
					"rounded" : 20
				}

			}
, 			{
				"box" : 				{
					"bgcolor" : [ 0.917647, 0.937255, 0.670588, 1.0 ],
					"id" : "obj-82",
					"maxclass" : "panel",
					"numinlets" : 1,
					"numoutlets" : 0,
					"patching_rect" : [ 435.0, 0.0, 428.0, 173.0 ],
					"presentation" : 1,
					"presentation_rect" : [ 5.0, 60.0, 385.0, 85.0 ],
					"rounded" : 20
				}

			}
 ],
		"lines" : [ 			{
				"patchline" : 				{
					"destination" : [ "obj-26", 0 ],
					"hidden" : 0,
					"midpoints" : [ 521.5, 519.0, 573.0, 519.0, 573.0, 606.0, 594.5, 606.0 ],
					"source" : [ "obj-53", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-25", 0 ],
					"hidden" : 0,
					"midpoints" : [ 521.5, 495.0, 594.5, 495.0 ],
					"source" : [ "obj-53", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-25", 2 ],
					"hidden" : 0,
					"midpoints" : [ 699.5, 495.0, 672.5, 495.0 ],
					"source" : [ "obj-55", 4 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-22", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-75", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-72", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-75", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-64", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-75", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-9", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-75", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-72", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-74", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-62", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1199.5, 315.0, 1134.5, 315.0 ],
					"source" : [ "obj-73", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-73", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-72", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-64", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-70", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-33", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1199.5, 252.0, 1257.0, 252.0, 1257.0, 315.0, 1249.5, 315.0 ],
					"source" : [ "obj-68", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-41", 1 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-62", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-29", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-60", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-21", 4 ],
					"hidden" : 0,
					"midpoints" : [ 739.5, 672.0, 793.5, 672.0 ],
					"source" : [ "obj-54", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-21", 3 ],
					"hidden" : 0,
					"midpoints" : [ 594.5, 672.0, 780.0, 672.0 ],
					"source" : [ "obj-26", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-21", 2 ],
					"hidden" : 0,
					"midpoints" : [ 594.5, 606.0, 571.0, 606.0, 571.0, 672.0, 766.5, 672.0 ],
					"source" : [ "obj-25", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-26", 2 ],
					"hidden" : 0,
					"midpoints" : [ 699.5, 570.0, 665.166687, 570.0 ],
					"source" : [ "obj-55", 4 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-26", 1 ],
					"hidden" : 0,
					"midpoints" : [ 673.25, 495.0, 693.0, 495.0, 693.0, 570.0, 629.833313, 570.0 ],
					"source" : [ "obj-55", 3 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-26", 3 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-27", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-26", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-25", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-25", 0 ],
					"hidden" : 0,
					"midpoints" : [ 594.5, 282.0, 594.5, 282.0 ],
					"source" : [ "obj-55", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-21", 1 ],
					"hidden" : 0,
					"midpoints" : [ 594.5, 246.0, 571.0, 246.0, 571.0, 672.0, 753.0, 672.0 ],
					"source" : [ "obj-29", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-21", 0 ],
					"hidden" : 0,
					"midpoints" : [ 929.5, 672.0, 739.5, 672.0 ],
					"source" : [ "obj-22", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-13", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-21", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-9", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-1", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-30", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1199.5, 99.0, 1176.0, 99.0, 1176.0, 117.0, 1071.0, 117.0, 1071.0, 285.0, 1229.5, 285.0 ],
					"source" : [ "obj-9", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-34", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-45", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-48", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1210.5, 491.0, 1175.5, 491.0 ],
					"source" : [ "obj-34", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-34", 1 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-33", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-45", 1 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-30", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-30", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-31", 1 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-40", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1106.5, 524.0, 1193.5, 524.0 ],
					"source" : [ "obj-41", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-42", 1 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-40", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-42", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-48", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-45", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-31", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-31", 0 ],
					"hidden" : 0,
					"midpoints" : [ 916.5, 323.0, 1209.5, 323.0 ],
					"source" : [ "obj-32", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-43", 1 ],
					"hidden" : 0,
					"midpoints" : [ 1158.5, 542.0, 1124.5, 542.0 ],
					"source" : [ "obj-38", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-38", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-44", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-56", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-43", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-56", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1175.5, 580.0, 1106.5, 580.0 ],
					"source" : [ "obj-42", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-43", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-41", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-41", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1038.5, 460.0, 1106.5, 460.0 ],
					"source" : [ "obj-39", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-35", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-37", 1 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-36", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-37", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-37", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1106.5, 527.0, 1045.5, 527.0 ],
					"source" : [ "obj-41", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-39", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-46", 1 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-32", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-46", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-55", 0 ],
					"hidden" : 0,
					"midpoints" : [ 521.5, 243.0, 594.5, 243.0 ],
					"source" : [ "obj-53", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-53", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-52", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-51", 2 ],
					"hidden" : 0,
					"midpoints" : [ 620.75, 348.0, 218.5, 348.0 ],
					"source" : [ "obj-55", 1 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-51", 1 ],
					"hidden" : 0,
					"midpoints" : [ 699.5, 348.0, 149.0, 348.0 ],
					"source" : [ "obj-55", 4 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-49", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-50", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-50", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-51", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-51", 0 ],
					"hidden" : 0,
					"midpoints" : [ 647.0, 348.0, 79.5, 348.0 ],
					"source" : [ "obj-55", 2 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-51", 0 ],
					"hidden" : 0,
					"midpoints" : [ 521.5, 348.0, 79.5, 348.0 ],
					"source" : [ "obj-53", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-46", 0 ],
					"hidden" : 0,
					"midpoints" : [ 594.5, 174.0, 916.5, 174.0 ],
					"source" : [ "obj-29", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-55", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-29", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-47", 0 ],
					"hidden" : 0,
					"midpoints" : [ 594.5, 156.0, 238.0, 156.0, 238.0, 417.0, 80.5, 417.0 ],
					"source" : [ "obj-29", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-2", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-29", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-29", 2 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-6", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-29", 1 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-8", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-8", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-10", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-8", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-12", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-15", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-13", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-15", 1 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-13", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-54", 0 ],
					"hidden" : 0,
					"midpoints" : [ 521.5, 426.0, 739.5, 426.0 ],
					"source" : [ "obj-53", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-54", 0 ],
					"hidden" : 0,
					"midpoints" : [ 1106.5, 616.0, 739.5, 616.0 ],
					"source" : [ "obj-56", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-54", 1 ],
					"hidden" : 0,
					"midpoints" : [ 673.25, 495.0, 774.833313, 495.0 ],
					"source" : [ "obj-55", 3 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-54", 2 ],
					"hidden" : 0,
					"midpoints" : [ 699.5, 570.0, 810.166687, 570.0 ],
					"source" : [ "obj-55", 4 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-54", 3 ],
					"hidden" : 0,
					"midpoints" : [ 620.75, 495.0, 845.5, 495.0 ],
					"source" : [ "obj-55", 1 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-52", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-3", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-10", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-7", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-12", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-11", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-3", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-5", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-3", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-75", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-67", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-64", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-16", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-67", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-24", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-16", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-68", 0 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-24", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-25", 1 ],
					"hidden" : 0,
					"midpoints" : [  ],
					"source" : [ "obj-58", 0 ]
				}

			}
, 			{
				"patchline" : 				{
					"destination" : [ "obj-58", 0 ],
					"hidden" : 0,
					"midpoints" : [ 647.0, 372.0, 634.5, 372.0 ],
					"source" : [ "obj-55", 2 ]
				}

			}
 ]
	}

}
