/*  Copyright 2018-2024 Lakshay Anand.
    All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    To view GNU General Public License  see <http://www.gnu.org/licenses/>.

    The chromoMap.js Javascript library depends on an open source software component.
    d3.js ,  https://github.com/d3/d3

  d3 license
----------------------------------------------------------------------
Copyright 2010-2024 Mike Bostock
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the author nor the names of contributors may be used to
  endorse or promote products derived from this software without specific prior
  written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------------
 */


/* chromoMap JavaScipt Library */

// importing savesvg and png 

// task 1: input files to json 

// task 2: input data transformation 

// task 3: annotating data logic 

// helper functions 

const print = (msg) => console.log(msg)

const validate = (para) => {
  para = typeof para !== "undefined" ? para : "default";
}

const ncol = (df) => Object.keys(df[0]).length
const nrow = (df) => df.length




/* Global variables */

var colors=[],plot_colors=[],chDataReduced,
    chr_spacing,y_val,ch_width,
    loci_width,top_margin,
    left_margin,x_scale_pos,
    chLinGradV=[],times,
    ch_curve,ch_gap,plot_spacing=[],
    v_align = false,plot_space,div_id,
    plot_shift_g = 0;

//var plot_y_labs=["#genes","-log10(FDR)"];
//var plot_y_labs=["#DEGs","logFC","Severe","Asymp.","logFC"];
//var plot_legend_label=["none","logFC","normalized expr.","normalized expr.","logFC"];
//var plot_legend_label=["logFC","normalized expr.","normalized expr."];
var make_plot_y_labs = true;

// creating the rounded ends of telomeric loci
function rounded_rect(x, y, w, h, r, tl, tr, bl, br) {
    var retval;
    retval  = "M" + (x + r) + "," + y;
    retval += "h" + (w - 2*r);
    if (tr) { retval += "a" + r + "," + r + " 0 0 1 " + r + "," + r; }
    else { retval += "h" + r; retval += "v" + r; }
    retval += "v" + (h - 2*r);
    if (br) { retval += "a" + r + "," + r + " 0 0 1 " + -r + "," + r; }
    else { retval += "v" + r; retval += "h" + -r; }
    retval += "h" + (2*r - w);
    if (bl) { retval += "a" + r + "," + r + " 0 0 1 " + -r + "," + -r; }
    else { retval += "h" + -r; retval += "v" + -r; }
    retval += "v" + (2*r - h);
    if (tl) { retval += "a" + r + "," + r + " 0 0 1 " + r + "," + -r; }
    else { retval += "v" + -r; retval += "h" + r; }
    retval += "z";
    return retval;
}

function getLongestChrLen(
  ch_files , ploidy = 1
){

  
  var ch_longest_len_arr = [];
  var longest_chr_info = [];
  var ch_data = [];

  for(let g = 0; g < ploidy; g++){
    ch_data.push(ch_files[g]);
    switch(ncol(ch_data[g])){
      case 3:
        ch_longest_len_arr.push(Math.max(...ch_data[g].map( (value) => value.end - value.start + 1 )));
        break;
      case 4:
        ch_longest_len_arr.push(Math.max(...ch_data[g].map( (value) => value.end - value.start + 1 )));
        cnt = true;
        break;
      default:
        print("Error: Allowed number of columns are 3 or 4.");

    }
  }

  return Math.max(...ch_longest_len_arr);
}



class chromoMapPlot {
  constructor(){
      
      
      
      
     
      // defaults properties 
      this.ploidy = 1;
      this.n_win_factor = 1;
      this.title = "";
      this.data_based_color_map = false;
      //this.data_type = ["numeric","categorical"];
      this.data_type = "numeric"
      this.plots_arr = ["none"];
      this.data_domain  = [];
      this.data_colors = [];
      this.heat_map = [true];
      this.legend = [true];
      return this;

      
  }

  

  getLongestLen(){
    
  }

  chromosomes(...ch_files){
    this.chrData = [...ch_files];
    this.chrDataCopy = JSON.parse(JSON.stringify(this.chrData));
    this.ch_longest_len = getLongestChrLen(this.chrData);
    //print(this.chrDataCopy)
    print(this.ch_longest_len);
    return this;
  }

  mappings(...data_files){
    this.annoData = [...data_files];
    this.annoDataCopy = JSON.parse(JSON.stringify(this.annoData));

    return this;
  }

  params(p){
    this.ploidy = p.ploidy;
    this.n_win_factor = p.n_win_factor;
    this.title = p.title;
    this.labels = p.labels;
    this.label_font = p.label_font;
    this.label_angle = p.label_angle;
    this.chr_curve = p.chr_curve;
    this.chr_length = p.chr_length;
    this.guides = p.guides;
    this.y_chr_scale = p.y_chr_scale;
    this.left_margin = p.left_margin;
    this.top_margin = p.top_margin;
    this.chr_text_x_shift = p.chr_text_x_shift;
    

    return this;

  }

  params_per_ploidy(p){
    this.data_colors = p.data_colors;
    this.anno_col = p.anno_col;
    this.data_domain = p.data_domain;
    this.display_chr = p.display_chr;
    this.chr_color = p.chr_color;
    this.chr_text = p.chr_text;

    return this;
  }

  createColorMapNumeric(){
    this.data_based_color_map = true;
    this.data_type = "numeric";
    return this;

  }

  colorMapCategorical(){
    this.data_based_color_map = true;
    this.data_type = "categorical";
    return this;

  }

  plots(...plots_arr){
    this.plots_arr = [...plots_arr];
    return this;
  }


  render(){
    $$$__chromoMap__$$$(this.chrData, this.annoData, {ploidy: this.ploidy,
                                                      n_win_factor : this.n_win_factor,
                                                      title : this.title,
                                                      data_based_color_map: this.data_based_color_map,
                                                      data_type: this.data_type,
                                                      plots: this.plots_arr,
                                                      data_domain: this.data_domain,
                                                      data_colors: this.data_colors,
                                                      heat_map: this.heat_map,
                                                      legend: this.legend  ,
                                                      labels: this.labels,
                                                      label_font: this.label_font,
                                                      label_angle: this.label_angle,
                                                      chr_curve: this.chr_curve,
                                                      chr_length: this.chr_length,
                                                      guides: this.guides,
                                                      anno_col: this.anno_col,
                                                      display_chr: this.display_chr,
                                                      chr_color: this.chr_color,
                                                      chr_text: this.chr_text,
                                                      y_chr_scale: this.y_chr_scale,
                                                      left_margin: this.left_margin,
                                                      top_margin: this.top_margin,
                                                      chr_text_x_shift: this.chr_text_x_shift
                                                    });
   }

  zoomView(factor = 0.03){

    this.chrData = JSON.parse(JSON.stringify(this.chrDataCopy));
    this.annoData = JSON.parse(JSON.stringify(this.annoDataCopy));
     // modify chdata
     for(let f = 0; f < this.ploidy; f++){
     this.chrData[f].forEach((x) => {
      if(x.start < this.ch_longest_len*factor) {
        x.start = this.ch_longest_len*factor;
      }; 
      if(x.end > (this.ch_longest_len - this.ch_longest_len*factor)){
        x.end = this.ch_longest_len - this.ch_longest_len*factor;
      } });}
     //print(this.annoData);
     this.render();
}

  set(){
   //set properties
}
}

function $$$__chromoMap__$$$(
                    /* Input files or objects */
                      // chromosome info 
                      ch_files = [],
                      // annotation info 
                      data_files = [],
                      //n_win_factor = 1,
                    /* multi-track or polyploidy */
                      //ploidy = 1,
                      {n_win_factor = 1 ,ploidy = 1,title = "",
                    /* general configuration  */
                      // for the plot and chromosomes 
                      
                      ch_gap = 5,
                      
                      top_margin = 30,
                      left_margin = 50,
                      chr_width = 20, //15
                      chr_length = 6, //4
                      chr_color = ["darkgray"],
                      canvas_width= null,
                      canvas_height = null,
                      anno_col = ["#3E9689"], //#10B85F
                      chr_text = true,
                      chr_text_x_shift = 0,
                      text_font_size = [15],
                      chr_curve = 5,
                      title_font_size = 15,
                      id = "chromap",

                      // scale 
                      chr_scale_ticks = 5,
                      y_chr_scale = 0,
                      scale_suffix = "bp",

                      // guides 
                      guides = false,
                      guides_color = "lightgrey",


                      interactivity = true,
                      display_chr = true,

                    /* data-based plots  */
                      // for heatmaps 
                      data_based_color_map ,
                      data_type ,

                      // data plot: scatter/bar/epi
                      plots = ["none"],

                      // add legends 
                      legend ,
                      // adjust legend 
                      lg_x = 0,
                      lg_y = 0,

                      // change heat colors 
                      data_colors = ["black"] ,
                      //discrete_domain = null,
                      //numeric_domain = null,
                      data_domain = [] ,
                      aggregate_func = ["avg"],
                      
                      tag_filter = [["none",0]],
                      plot_height =[30],
                      plot_ticks = [4],
                      plot_color = ["blue"],
                      plot_y_domain = [[0,0]],
                      ref_line = [false],
                      refl_pos = [0],
                      refl_color = ["grey"],
                      refl_stroke_w = [2],
                      tagColor = ["red"],
                      heat_map ,
                      plot_shift = [1],
                      plot_legend_labels = [],
                      cat_legend_label = [],
                      plot_y_labels = [],
                      plot_y_lab_x = 10,
                      plot_y_lab_y = 0,
                      plot_y_lab_size = 15,

                      // filters
                      plot_filter = [["none",0]],

                    /* segmental annotaion */   
                      segment_annotation = false,
                    /* labelling */
                      labels = false,
                      label_font = 9,
                      label_angle = -90,
                    /* Grid lines highliting */
                      
                      vertical_grid = false,
                      grid_array = [0,5000,10000],
                      grid_color = "grey",
                      grid_text = null,
                      grid_text_size = 12,
                      grid_text_y = 20,
                    
                    /* hyperlinking */
                      hlinks = false,
                    /* quick zoom-in */
                     region = null,

                    /* high resolution */
                      //n_win_factor = 1,
                      fixed_window = false,
                      window_size = null,
                      win_summary_display = false,
                      remove_last_window = true,

                    /* 2D-chromosome plots */
                      chr_2D_plot = false,
                      ch2D_colors = null,
                      ch2D_cat_order = null,
                      ch2D_lg_x = 0,
                      ch2D_lg_y = 0,
                      ann_h = 1,

                    /* chromLinks */
                      //data

                      show_links = false,
                      loci_links = "none",
                      directed_edges = false,
                      
                      links_colors = null,
                      links_lg_x = 0,
                      links_lg_y = 0,

                    /* exporting options */
                      export_options = true,
                    } = {}

                   
                   
                    
                     
                      
) {
   // chromoMap() START

   
   print(data_files);

   if(display_chr === true){
    display_chr = Array(ploidy).fill(true)
   }

   if(chr_text === true){
    chr_text = Array(ploidy).fill(true)
   }

   /* Input Validation */ 

  var color_map = data_based_color_map;
  //var color_scale = data_type;
  //color_scale =color_scale[0];


  

  /* feching the data to render chromoMap */
  const chr_data = [];
  const max_ch_domain = []; 
  var ch_longest_len_arr = [];
  var longest_chr_info = [];
  const ch_data = [];
  var cnt = false;


  for(let g = 0; g < ploidy; g++){
    ch_data.push(ch_files[g]);
    switch(ncol(ch_data[g])){
      case 3:
        ch_longest_len_arr.push(Math.max(...ch_data[g].map( (value) => value.end - value.start + 1 )));
        break;
      case 4:
        ch_longest_len_arr.push(Math.max(...ch_data[g].map( (value) => value.end - value.start + 1 )));
        cnt = true;
        break;
      default:
        print("Error: Allowed number of columns are 3 or 4.");

    }
  }

  var ch_longest_len = Math.max(...ch_longest_len_arr);

  //if(getLongest) return ch_longest_len;

  /*
  function zoom_update(factor) {
    data.forEach((x) => {if(x.start < ch_longest_len*factor) {x.start = ch_longest_len*factor;}; if(x.stop > ch_longest_len - ch_longest_len*factor){x.stop = ch_longest_len - ch_longest_len*factor}});

  }  */
  

  for(let g = 0; g < ploidy; g++){
         // for start
    chr_data.push(ch_files[g]);
       
    print(`Number of Chromosomes in set ${g}: ${nrow(chr_data[g])}\n`);
    
    if(!fixed_window){

      switch(ncol(chr_data[g])){
        case 3:
          chr_data[g].forEach((value) => {
            //value.real_start = value.start;
            //value.start = 1;
            let sq_len = value.end - value.start +1;
            let sq_len2 = value.start - 1;
            value.ch_average = (sq_len)/ch_longest_len;
            value.n = (((sq_len)/ch_longest_len).toFixed(2))*(100*n_win_factor);
            value.n0 = (((sq_len2)/ch_longest_len).toFixed(2))*(100*n_win_factor);
            //value.n = value.n2 + value.n1;
            //value.n0 = 0;
            value.seq_len = sq_len;
            value.ploidy = g;
          });
          break;
        case 4:
          chr_data[g].forEach((value) => {
            let sq_len = value.end - value.start +1;
            value.ch_average = (sq_len)/ch_longest_len;
            value.n = (((sq_len)/ch_longest_len).toFixed(2))*(100*n_win_factor);
            value.seq_len = sq_len;
            value.cnt_proprtion = (value.cnt_start/(sq_len)).toFixed(2);
            value.p = (((value.cnt_start/sq_len).toFixed(2))*((((sq_len)/ch_longest_len).toFixed(2))*(100*n_win_factor))).toFixed();
            value.q = ((((sq_len)/ch_longest_len).toFixed(2))*(100*n_win_factor)).toFixed() - (((value.cnt_start/sq_len).toFixed(2))*((((sq_len)/ch_longest_len).toFixed(2))*(100*n_win_factor))).toFixed();
            value.ploidy = g;
          });
          break;
        default:
          print("Error: Allowed number of columns are 3 or 4.");
      }

    } else {

      switch(ncol(chr_data[g])){
        case 3:
          chr_data[g].forEach((value) => {
            let sq_len = value.end - value.start +1;
            value.ch_average = (sq_len)/ch_longest_len;
            value.n = Math.ceil(seq_len/window_size);
            value.seq_len = sq_len;
            value.ploidy = g;
          });
          break;
        case 4:
          chr_data[g].forEach((value) => {
            let sq_len = value.end - value.start +1;
            value.ch_average = (sq_len)/ch_longest_len;
            value.n = Math.ceil(seq_len/window_size);
            value.seq_len = sq_len;
            value.cnt_proprtion = (value.cnt_start/(sq_len)).toFixed(2);
            value.p = Math.ceil(value.cnt_start/window_size);
            value.q = Math.ceil((sq_len - value.cnt_start)/window_size);
            value.ploidy = g;
          });
          break;
        default:
          print("Error: Allowed number of columns are 3 or 4.");
      }

      
    }
    
    longest_chr_info.push(chr_data[g].filter((value) => value.seq_len == ch_longest_len_arr[g])[0]);

    
   //for end
  }

  var long_chr = longest_chr_info.filter((value) => value.seq_len == ch_longest_len);

  /* Creating genomic ranges */

  var mega_list_of_ranges = [];
  var ch_name_list = [];
  var tempInput = [];

  //check if cnt
  if(fixed_window){
    if(cnt){
      print("Fixed window visualization is not supported with centromeres.");
    }
  }

  for(let g = 0; g < ploidy; g++){
  
    var list_of_ranges = {};
    var namee = [];
    let nr = nrow(chr_data[g]);

    chr_data[g].forEach((value) => {
      var ch_start = value.start;
      var ch_end = value.end;
      var ch_loci = value.n;

      if(!fixed_window){       
        var step_size = (ch_end - ch_start + 1)/ch_loci;
      } else {
        var step_size = window_size;
      }

      var shift_distance = value.start - 1;

      var brkpoints = Array.from({length: ch_loci}, (_, i) => i + 1).map((x) => Math.ceil(x*step_size));

      if(fixed_window){
        if(remove_last_window){
          if(((ch_end-ch_start+1) - brkpoints[(brkpoints.length - 2)]) != step_size){
            
            brkpoints = brkpoints.pop();
            value.n = value.n - 1;
            value.end = brkpoints[brkpoints.length];
            
          }}
      }

      var range_start = [];
      var range_end = [];
      var ch_name = [];
      var ch_loci_end;

      for(let i = 0; i < brkpoints.length; i++){
        
        
        var ch_loci_start = ch_start;
        if(i != brkpoints.length){
          ch_loci_end = brkpoints[i] + shift_distance;
        } else {
          if(!remove_last_window){
            ch_loci_end = ch_end;
          } else {
            ch_loci_end = brkpoints[i] + shift_distance;
          }
        }
        ch_start = ch_loci_end + 1;
        
        
        range_start.push(ch_loci_start);
        range_end.push(ch_loci_end);
        ch_name.push(value.name);
        
      }

      var range_df = [];
      
      for(let r = 0; r < range_start.length; r++){
        range_df.push( {
          "range_start" : range_start[r],
          "range_end" : range_end[r],
          "ch_name" : ch_name[r],
          "idx" : r + 1
        });
      }

      ch_name = ch_name.filter((x, i, a) => a.indexOf(x) === i);
      list_of_ranges[value.name] = range_df;
      namee.push(ch_name);


    });

    
    ch_name_list.push(namee);
    //names(list_of_ranges)=chr_data[[g]]$ch_name
    mega_list_of_ranges.push(list_of_ranges);
    
  }

  /* Window summary display */


  /* calculating ch domain */

  var ch_domain = [];
  //#print(mega.list.of.ranges[[as.numeric(long.chr[1,2])]])
  
  var tmp_list_obj = mega_list_of_ranges[long_chr[0].ploidy];
  var longest_chr_df = tmp_list_obj[long_chr[0].name]

  //print(nrow(longest_chr_df));
  

  // adjusting canvas 

  var dom_idx = 0;
  for(let f = 0; f < nrow(longest_chr_df);f++){
    if(f != (nrow(longest_chr_df)-1)){
      ch_domain.push(longest_chr_df[f].range_start);
      dom_idx = dom_idx + 1;
    } else {
      ch_domain.push(longest_chr_df[f].range_start);
      dom_idx = dom_idx + 1;
      ch_domain.push(longest_chr_df[f].range_end);
      
    }

  }
  
  
  //ch.domain <- as.character(ch.domain)
  //#print(ch.domain)

  //adjusting for width
  print(long_chr);
  if(canvas_width == null){
    if(!fixed_window){
      canvas_width = 100*n_win_factor*chr_length + left_margin + 100
    } else {
      canvas_width = chr_length*nrow(longest_chr_df) + left_margin + 100
    }
  }
  
  //defining the x-scale
  var win_scale;
  if(!fixed_window){
    win_scale = 100*n_win_factor;
  } else {
    win_scale = nrow(longest_chr_df);
  }


  /* Assigning Loci/geneomic ranges */
   print(chr_data);
  //assigning loc for each gene or elemnt
  print("Processing data.. \n")
  
  var inputData = [];
  var labels_ids = [];
  
  if(!segment_annotation){
  
  for(let h = 0; h < ploidy; h++){


    inputData.push(data_files[h]);

    var temp_outer = [], temp_dff;
    var bf_rw = nrow(inputData[h]);

    chr_data[h].forEach((value) => {
        // check if chr exist 
        temp_dff = inputData[h].filter((d) => d.ch_name == value.name);
        temp_dff = temp_dff.filter((c) => c.ch_start >= value.start && c.ch_end <= value.end);
        temp_outer = [...temp_outer,...temp_dff];
    });

   
    inputData[h] = temp_outer;
    var af_rw = nrow(inputData[h]);
    //var data_col2 = ncol(inputData[h]);
    if(bf_rw != af_rw){
      print(`WARNING: $(bf_rw-af_rw) out-of-bound annotations are removed in chromosome set $(h) \n`)
    }

    // for 2D plots 
    if(chr_2D_plot[h]){
      color_map = true;
      plots[h] = "2d";
      data_based_color_map = true;
      color_scale = "linear";
      heat_map[h] = false;
      legend[h] = false;

    }

    //colnames(inputData[[h]])=c("name","ch_name","ch_start","ch_end","data","hlink","cate")
    
    


    //print("Number of annotations in data set ",h,":",nrow(inputData[[h]]),"\n");

    var loci = [];
    var chr_id = []
    var loci_start = [];
    var loci_end = [];
    var label = [];
    var chrom = [];

    inputData[h].forEach((value) => {
      var temp_list = mega_list_of_ranges[h];
      //names(temp.list)=ch.name.list[[h]]
      
      var temp_df = temp_list[value.ch_name];

      temp_df = temp_df.filter((d) => d.range_start <= value.ch_start && d.range_end >= value.ch_start);
      value.loci = id+"-"+value.ch_name+"-"+temp_df[0].idx+"-"+(h+1);
      value.chr_id = id+"-"+value.ch_name+"-"+(h+1);
      value.loci_start = temp_df[0].range_start;
      value.loci_end = temp_df[0].range_end;
      value.label = id+"-L"+value.ch_name+"-"+temp_df[0].idx+"-"+(h+1);
      value.chrom = value.ch_name;


    });
    
  }
  
  } else {
    // segmental annotations 
  }



  // Data domains for heatmap 

  // compute data domains from data if data_domain is null.

 

  if(color_map){

    if(data_domain.length == 0){

    if( data_type == "numeric"){

      for(let h = 0; h < ploidy; h++){

        let dd = d3.extent(inputData[h], (d) => { return d.data;});

        if(dd[0] < 0 && dd[1] > 0){

          data_domain.push([dd[0],0,dd[1]]);

        } else {

        data_domain.push(dd);
                }

      }


    } 

    if( data_type == "categorical"){


    }
} 
  }


  // assigning data colors 
 // selecting schemes 

 

 console.log(data_colors.length)

 if(color_map){

  for(let h = 0; h < ploidy; h++){

  
  if(data_colors.length == 0) data_colors[h] = d3.schemeAccent;

  data_colors[h] = data_colors[h].slice(0,data_domain[h].length)

}

 }


   
  print(inputData);
  print("this is where ")
  print(data_colors)
  print(data_domain)


  $_chromoMap$(inputData,chr_data,
    ploidy_n=ploidy,
    title=title,
    cnt=cnt,
    ch_gap=ch_gap,
    top_margin=top_margin,
    left_margin=left_margin,
    chr_width=chr_width,
    chr_length=chr_length,
    chr_col=chr_color,
    heatmap=color_map,
    ch_domain=ch_domain,
    lg_x=lg_x,
    lg_y=lg_y,
    heat_scale=data_type,
    labels=labels,
    div_id=id,
    w=canvas_width,
    h=canvas_height,
    rng=data_domain,
    heat_col=data_colors,
    an_col=anno_col,
    ch_text=chr_text,
    legend=legend,
    aggregate_func=aggregate_func,
    plots=plots,
    tag_filter = tag_filter,
    plot_height = plot_height,
    plot_ticks = plot_ticks,
    plot_color = plot_color,
    plot_y_domain = plot_y_domain,
    ref_line = ref_line,
    refl_pos = refl_pos,
    refl_color = refl_color,
    refl_stroke_w = refl_stroke_w,
    tagColor = tagColor,
    renderHeat = heat_map,
    text_font_size = text_font_size,
    chr_curve = chr_curve,
    title_font_size = title_font_size,
    label_font = label_font,
    label_angle = label_angle,
    vertical_grid = vertical_grid,
    grid_array = grid_array,
    grid_color = grid_color,
    plot_filter = plot_filter,
    loci_links = [],
    uniq_cates = [],
    scatter_col = ch2D_colors,
    grid_text = grid_text,
    grid_text_size = grid_text_size,
    grid_text_y = grid_text_y,
    scatter_mapping = false,
    scatter_lg_x = ch2D_lg_x,
    scatter_lg_y = ch2D_lg_y,
    show_links = show_links,
    seg_anno = segment_annotation,
    directed_edges = directed_edges,
    y_chr_scale = y_chr_scale,
    links_colors = links_colors,
    links_lg_x = links_lg_x,
    links_lg_y = links_lg_y,
    links_color_maps=false,
    win_scale = win_scale,
    scale_ticks = chr_scale_ticks,
    export_options = export_options,
    guides = guides,
    guides_color = guides_color,
    ann_h = ann_h,
    display_chr = display_chr,
    plot_shift = plot_shift,
    plot_legend_label = plot_legend_labels,
    cat_legend_lab = cat_legend_label,
    plot_y_labs = plot_y_labels,
    plot_y_lab_x = plot_y_lab_x,
    plot_y_lab_y = plot_y_lab_y,
    plot_y_lab_size = plot_y_lab_size,
    scale_suffix = scale_suffix,
    interactivity = interactivity,
    chr_text_x_shift = chr_text_x_shift);

  // chromoMap() End 
}


//LOGIC FOR THE CHROMOMAP
function $_chromoMap$(chData,nLoci,ploidy_n,title,cnt,ch_gap,
  top_margin,left_margin,chr_width,chr_length,chr_col,heatmap,
  ch_domain,lg_x,lg_y,heat_scale,
  labels,div_id,width,height,rng,
  heat_col,an_col,ch_text,legend,aggregate_func,plots,tag_filter,
  plot_h,plot_ticks,plot_color,plot_y_domain,
  ref_line,refl_pos,refl_color,refl_stroke_w,
  tagColor,renderHeat,text_font_size,chr_curve,title_font_size,
  label_font,label_angle, grid_array,vertical_grid,grid_color,
  plot_filter,loci_links,uniq_cates,scatter_col,
  grid_text,grid_text_size,grid_text_y,
  scatter_mapping,scatter_lg_x,scatter_lg_y,
  show_links,seg_anno,directed_edges,y_chr_scale,
  links_colors,links_lg_x,links_lg_y,links_color_maps,
  win_scale,scale_ticks,export_options,
  guides,guides_color,ann_h,display_chr,plot_shift,
  plot_legend_label,cat_legend_lab,plot_y_labs,
  plot_y_lab_x,plot_y_lab_y,plot_y_lab_size,
  scale_suffix,interactivity,chr_text_x_shift){




//swapping margins
if(v_align){
var t=top_margin;
top_margin=left_margin;
left_margin=t;

}




if(!labels){
y_val = top_margin;
} else {

y_val = top_margin+10;
}




ch_width = chr_width; /*height of chromosome bar */
arc_radius=chr_width/2;

loci_width=chr_length; /*widht of each loci determined for ch length */

ch_curve=chr_curve; /* the curve at the end loci*/

var plot_height = plot_h;
//console.log(plot_h);
if(plots[0] == "none"){
var plot_spacing = [];
for(var h=0;h<ploidy_n;h++){
plot_spacing.push(0);}
plot_space = 0

} else {

plot_spacing = plot_height.map( function(value) {
return value + 18; } );


plot_space = plot_spacing.reduce(function(a, b){
return a + b;
}, 0);

//plot_space = plot_space + (ploidy_n - 1)*(ch_gap*3);
}
//console.log(plot_space);
if(labels){
plot_spacing = plot_spacing.map( function(value) {
return value + 50; } );

plot_space = plot_spacing.reduce(function(a, b){
return a + b;
}, 0);
}

var plot_padding = 15;

if(!labels){
chr_spacing= ploidy_n*chr_width + ploidy_n + ch_gap*2 + plot_space;
} else {
chr_spacing= ploidy_n*chr_width + ploidy_n + ch_gap*2 + plot_space;
}



ttl=title;


//     var arc = d3.arc()
//  .innerRadius(0)
//  .outerRadius(arc_radius)
//  .startAngle(0)
//  .endAngle(Math.PI);

//  var arcLeft = d3.arc()
//  .innerRadius(0)
//  .outerRadius(arc_radius)
//  .startAngle(Math.PI)
//  .endAngle(2*Math.PI);

//  var arcTop = d3.arc()
//  .innerRadius(0)
//  .outerRadius(arc_radius)
//  .startAngle(-Math.PI/2)
//  .endAngle(Math.PI/2);

//  var arcBot = d3.arc()
//  .innerRadius(0)
//  .outerRadius(arc_radius)
//  .startAngle(Math.PI/2)
//  .endAngle((3*Math.PI)/2);

// automaticlally update height if not given

if(height == null){
if(!legend.includes(true)){
height = ( top_margin+(chr_spacing)*nLoci[0].length+10)- y_chr_scale + 100;
} else {
height = ( top_margin+(chr_spacing)*nLoci[0].length+10)- y_chr_scale + 300;
}
}
/*
let cmap_divvv = document.getElementById(div_id);
if(cmap_divvv !== null) {
  let f1 = document.getElementById('cmap_canvas');
  f1.removeChild(f1.lastElementChild);
  let f2 = document.getElementById('dwopt');
  if(f2 !== null) f2.removeChild(f2.lastElementChild);
  let f3 = document.getElementById(div_id+"-svg");
  if(f3 !== null) f3.removeChild(f3.lastElementChild);
}  */



 

// adding download buttons
if(export_options){
var css = `
#save_png2 {
background-color: #4CAF50; 
color: white;
padding: 15px 32px;
text-align: center;
font-size: 16px;
width:100%;
border:solid;
border-radius:10px;
}
#save_png2:hover{
background-color: white;
color: black;
cursor:pointer;
}
#save_svg2 {
background-color: #4CAF50; 
color: white;
padding: 15px 32px;
text-align: center;
font-size: 16px;
width:100%;
border:solid;
border-radius:10px;
}
#save_svg2:hover{
background-color: white;
color: black;
cursor:pointer;
}

#dwopt {
  background: red;
  display: flex;
}


.chromoMap_toolbar {
  width: 100%;
  height: 130px;
  background: transparent;
  display: flex;
  flex-direction: column;
  border-radius: 10px;
  justify-content: center;
  align-items: center;
}
.navi{
   position: relative;
   width: 100%;
   height: 30px;
   background: #C0C0C0;
   display: flex;
   justify-content: center;
   align-items: center;
   border-top-left-radius: 10px;
   border-top-right-radius: 10px;
}

.navi ul {
   display: flex;
   width: 100%;
   padding: 0;
   
   
}

.navi ul li {
   list-style: none;
   position: relative;
   flex-grow: 1;
   height: 30px;
   z-index: 2;
   transition: 0.2s;
   border-top-left-radius: 10px;
   border-top-right-radius: 10px;
  
   
}



.navi ul li a {
   display: flex;
   justify-content: center;
   align-items: center;
   text-decoration: none;
   color: #36454F;
}

.navi ul li a .text {
   
  
   color: #36454F;
   padding: 2px 7px;
   
   font-weight: 400;
   font-size: 0.75em;
   letter-spacing: 0.2em;
   transition: 0.5s;
   transform: translateY(5px);
   

}
.navi ul li.active a .text {
   color: #C0C0C0;
}

.navi ul li.active {
   background: #36454F;
}

.indicator{
   position: absolute;
   width: 60px;
   height: 60px;
   background: #36454F;
   border-radius: 50%;
   top:-35px;
   left:10vw;
   z-index: 1;
   transition: 0.2s;

}

.indicator::before{
   content: '';
   position: absolute;
   width: 30px;
   height: 30px;
   background: transparent;
   border-radius: 50%;
   top:5px;
   left:-28px;
   box-shadow: 15px 18px #36454F;
}

.indicator::after{
   content: '';
   position: absolute;
   width: 30px;
   height: 30px;
   background: transparent;
   border-radius: 50%;
   top:5px;
   right:-28px;
   box-shadow: -15px 18px #36454F;
}

.navi ul li:nth-child(1).active ~ .indicator {
   transform: translateX(calc(25vw * 0));
}

.navi ul li:nth-child(2).active ~ .indicator {
   transform: translateX(calc(25vw * 1));
}

.navi ul li:nth-child(3).active ~ .indicator {
   transform: translateX(calc(25vw * 2));
}
.navi ul li:nth-child(4).active ~ .indicator {
   transform: translateX(calc(25vw * 3));
}



.chromoMap_toolarea {
   background: #E5E4E2;
   height: 70px;
   width: 100%;
   
   border-bottom-left-radius: 10px;
   border-bottom-right-radius: 10px;
}

[data-tab-content]{
   display: none;
   transition: 0.5s;
}

.visible[data-tab-content]{
   display: block;
}

.chromoMap_toolarea .panel {
   width: 100%;
   height: 100%;
   display: flex;
   justify-content: space-evenly;
   align-items: center;
   padding-top: 5px;
   gap:5px;
   transition: 0.5s;
}

.chromoMap_toolarea .panel .divider {
   width: 2px;
   height: 60px;
   position: relative;
   
   background: #D3D3D3;
   box-shadow: 1px 0px 1px #D3D3D3;
}

.chromoMap_toolarea .panel .area {
   display: flex;
   flex-direction: row;
   justify-content: flex-start;
   
}



.chromoMap_toolarea .panel .area p{
   transform: rotate(-90deg);
   color: #36454F;
   padding: 2px 7px;
   
   font-weight: 400;
   font-size: 0.75em;
   letter-spacing: 0.05em;
   transition: 0.5s;
   word-break: break-all;
   white-space: normal;
}

#chr_SVG{
   stroke: url(#lg1);
   fill: url(#lg1);
   stroke-width: 0;
   height: 60px;
   width: 60px;
   transform: scale(1);
  
   
   

  
            
}

#chr_1 {
   transform: translateX(45px);
   
}

`;

var ch_tool_bar = `

<div style="width:100%; height:100%" id="cmap_canvas">
</div>
    <div class="chromoMap_toolbar">

    
        <div class="navi">
            <ul>
                <li class="list active" data-tab-target="#chromoMap_panel001">
                    <a href="#">
                        <span class="text"> Navigation </span>
                    </a>
                </li>
                <li class="list" data-tab-target="#chromoMap_panel002">
                    <a href="#">
                        <span class="text"> Properties </span>
                    </a>
                </li>
                <li class="list" data-tab-target="#chromoMap_panel003">
                    <a href="#">
                        <span class="text"> Markers </span>
                    </a>
                </li>
                <li class="list" data-tab-target="#chromoMap_panel004">
                    <a href="#">
                        <span class="text"> Export </span>
                    </a>
                </li>
                <div class="indicator">
                    <svg id="chr_SVG" viewBox="0 0 600 480" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg">
                        <!-- Created with SVG-edit - http://svg-edit.googlecode.com/ -->
                        <defs>
                            <linearGradient id="lg1">
                                <stop offset="0%" stop-color="#d2305c"></stop>
                                <stop offset="25%" stop-color="#ff5e00"></stop>
                                <stop offset="50%" stop-color="#ffeb00"></stop>
                                <stop offset="75%" stop-color="#1cff00"></stop>
                                <stop offset="100%" stop-color="#00b9ff"></stop>
                            </linearGradient>
                        </defs>
                        <g>
                         
                            <path id="chr_1" d="m94,164c-0.7045,-0.7045 -1.3069,-1.29317 -2,-2c-0.52612,-0.53654 -1.17366,-1.32834 -1.65717,-2.01433c-0.36905,-0.52359 -0.94214,-1.38969 -1.543,-2.28568c-0.40298,-0.60091 -0.8091,-1.20859 -1.257,-1.8c-0.69189,-0.9136 -1.23023,-1.50334 -1.7895,-2.13199c-0.60078,-0.67531 -1.25092,-1.41533 -1.92,-2.26801c-0.36061,-0.45956 -1.08768,-1.4613 -1.45867,-1.98666c-0.37637,-0.53297 -0.74805,-1.07422 -1.49066,-2.14667c-0.71611,-1.03419 -1.05985,-1.52528 -1.72234,-2.44783c-0.31383,-0.43703 -0.9157,-1.26099 -1.2,-1.65268c-0.79698,-1.09805 -1.27895,-1.77902 -1.72567,-2.43016c-0.42784,-0.62363 -1.23325,-1.83784 -1.63734,-2.436c-0.40658,-0.60185 -0.82199,-1.2097 -1.284,-1.79999c-0.71667,-0.91565 -1.27222,-1.50018 -1.81467,-2.12083c-0.55235,-0.63199 -1.34888,-1.63565 -1.816,-2.3645c-0.47468,-0.74065 -0.67829,-1.13882 -1.13667,-1.91917c-0.70526,-1.20067 -1.23077,-1.99721 -1.5,-2.39999c-0.54618,-0.81711 -0.84289,-1.21854 -1.36333,-2.08083c-0.79488,-1.31699 -1.00799,-1.79427 -1.448,-2.74c-0.21985,-0.47252 -0.62994,-1.42268 -0.83467,-1.89067c-0.59697,-1.36465 -1.00282,-2.2328 -1.212,-2.65867c-0.41931,-0.85365 -0.86915,-1.71236 -1.356,-2.592c-0.51023,-0.92188 -0.78944,-1.38617 -1.3125,-2.35417c-0.26352,-0.48768 -0.55668,-0.96383 -1.1875,-2.47916c-0.20557,-0.49382 -0.53002,-1.49361 -0.8125,-2.47916c-0.13871,-0.48392 -0.38216,-1.44238 -0.6875,-2.35417c-0.298,-0.88985 -0.66468,-1.7434 -1.044,-2.592c-0.19149,-0.42839 -0.79501,-1.73889 -1.188,-2.65867c-0.19973,-0.46745 -0.584,-1.41465 -0.76534,-1.89067c-0.35759,-0.93872 -0.81722,-2.3032 -0.952,-2.74133c-0.26367,-0.85714 -0.5019,-1.7116 -0.768,-2.56799c-0.27501,-0.88509 -0.59777,-1.80081 -0.94133,-2.74133c-0.17354,-0.47508 -0.52631,-1.42048 -0.69067,-1.89067c-0.47482,-1.35836 -0.74134,-2.22665 -0.94233,-3.084c-0.19891,-0.84846 -0.28809,-1.28032 -0.42033,-2.62133c-0.09165,-0.92935 -0.11922,-1.87878 -0.14667,-2.35334c-0.05447,-0.94158 -0.18734,-2.30519 -0.32733,-3.16816c-0.06844,-0.42188 -0.23752,-1.24395 -0.36333,-2.05584c-0.18563,-1.19783 -0.23065,-1.99688 -0.284,-2.78049c-0.052,-0.76372 -0.20966,-2.20175 -0.35266,-2.86367c-0.2053,-0.95026 -0.41809,-1.54804 -0.562,-2.75583c-0.10446,-0.87666 -0.08308,-1.44416 -0.08533,-2.71467c-0.00114,-0.645 -0.17862,-1.70415 0.08533,-2.48534c0.30055,-0.88952 0.70475,-1.84552 1.4355,-2.1805c0.79971,-0.36659 1.73501,-0.46643 2.3645,-0.73683c0.91221,-0.39185 1.80133,-0.28333 2.8,-0.284c0.8355,-0.00056 1.41917,0 2.6145,0c0.8895,0 1.743,0 2.77883,0c0.696,0 1.60542,-0.01155 2.604,0.00134c0.79667,0.01028 1.59821,-0.05371 2.5075,0.35583c0.75013,0.33786 1.47598,0.92028 2.15717,1.34283c0.91917,0.57018 1.52338,1.12104 2.30533,1.6c0.68356,0.4187 1.60711,0.98178 2.146,1.5c0.72492,0.69711 1.29868,1.52291 1.70133,2.2c0.49673,0.83529 1.09586,1.45288 1.37867,2.32133c0.21347,0.65554 0.32869,1.66634 0.82583,2.43134c0.48196,0.74165 0.90187,1.4548 1.2955,2.24733c0.39577,0.79684 0.75948,1.59058 0.94733,2.5955c0.13869,0.74187 0.4433,1.63667 0.854,2.37383c0.43503,0.78083 0.69024,1.55473 0.99867,2.34801c0.3894,1.00152 0.70039,1.66164 1.2,2.68266c0.29574,0.60439 0.70158,1.44467 0.87867,2.42133c0.14446,0.79672 0.10179,1.64045 0.34283,2.65717c0.16569,0.69891 0.58549,1.51738 0.9785,2.32284c0.40587,0.83183 0.80563,1.76321 1.10001,2.55133c0.2039,0.54589 0.4831,1.37026 0.89999,2.36333c0.38453,0.91599 0.80862,1.84744 1.10001,2.64117c0.20612,0.56145 0.50173,1.4436 0.8,2.34299c0.30176,0.90993 0.51257,1.53212 0.9,2.8665c0.2081,0.71674 0.40201,1.47191 0.636,2.244c0.24046,0.79341 0.68127,1.99126 0.84933,2.39066c0.34209,0.813 0.72801,1.63344 1.116,2.484c0.40144,0.88007 0.79946,1.79644 0.99416,2.26283c0.38792,0.92924 0.74958,1.84113 0.91917,2.2825c0.47023,1.22379 0.70786,1.96062 1.11067,2.912c0.24381,0.57584 0.84737,1.58092 1.37466,2.2365c0.50532,0.62824 1.13939,1.40668 1.58933,2.22083c0.43497,0.78706 0.81068,1.60001 1.21067,2.4c0.4,0.8 0.8,1.6 1.2,2.4c0.4,0.8 0.8,1.60001 1.2,2.39999c0.4,0.8 0.8,1.60001 1.2,2.40001c0.3,0.59999 0.8,1.59999 1.2,2.39999c0.4,0.8 0.8,1.60001 1.1,2.2c0.5,1 0.79731,1.60136 1.3,2.60001c0.39787,0.79041 0.72932,1.56909 1.3045,2.24283c0.52199,0.61143 1.36263,1.03304 2.04284,1.45717c0.91151,0.56836 1.12923,1.55948 1.81933,2.13333c0.65444,0.5442 1.74945,0.07587 2.22267,-0.72665c0.34493,-0.58496 0.49967,-1.74887 0.558,-2.42551c0.06403,-0.74278 0.07961,-1.54985 0.21933,-2.84784c0.09912,-0.92079 0.18747,-1.39026 0.33334,-2.35416c0.14883,-0.98344 0.24286,-1.47783 0.3785,-2.97917c0.0451,-0.49927 0.092,-1.49986 0.11083,-2.5c0.00941,-0.49995 0.01067,-2 0.01067,-2.50017c0,-1.00433 0,-1.5105 0,-2.53584c0,-1.04933 0,-2.13066 0,-2.68533c0,-1.13133 0,-1.70399 0,-2.85333c0,-0.57266 0.00172,-1.14133 -0.00017,-2.804c-0.00061,-0.53733 -0.00652,-1.59609 -0.02067,-2.12133c-0.02836,-1.05291 -0.04014,-1.5862 -0.20066,-3.2215c-0.05544,-0.56472 -0.21243,-1.71481 -0.29934,-2.29934c-0.08774,-0.59014 -0.27993,-1.78059 -0.47916,-2.97916c-0.09972,-0.59988 -0.29967,-1.79556 -0.4,-2.38934c-0.09967,-0.58988 -0.29891,-1.75368 -0.4,-2.32533c-0.19788,-1.11904 -0.39954,-2.20142 -0.5,-2.73267c-0.1991,-1.05283 -0.30214,-1.57759 -0.5,-2.63667c-0.20219,-1.08226 -0.40058,-2.20123 -0.5,-2.77317c-0.10059,-0.57872 -0.30024,-1.75346 -0.39999,-2.34733c-0.10025,-0.59679 -0.3992,-2.39567 -0.50134,-2.9955c-0.1022,-0.60017 -0.31662,-1.80057 -0.43467,-2.4c-0.11827,-0.60057 -0.52727,-2.40089 -0.68517,-2.99983c-0.1581,-0.59975 -0.3308,-1.19461 -0.658,-2.37933c-0.16138,-0.58432 -0.3576,-1.15396 -0.6875,-2.85417c-0.10529,-0.54265 -0.22742,-1.61339 -0.26134,-2.14133c-0.06743,-1.04971 -0.0675,-1.57798 -0.072,-3.19199c-0.00155,-0.55483 0,-1.69051 0,-2.26933c0,-0.58484 0,-1.7685 0,-2.36533c0,-1.19867 0,-2.39867 0,-2.99867c0,-0.6 0,-1.8 0,-2.4c0,-0.6 0,-2.39983 0,-2.99866c0,-0.59684 0,-1.7805 0,-2.36533c0,-1.15067 0,-2.26933 0,-3.34933c0,-0.525 0,-1.545 0,-2.04c0,-0.96 0,-2.32516 0,-3.1895c0,-0.42117 0,-1.2465 0,-2.05583c0,-1.20134 0,-2.00134 0,-2.80133c0,-0.8 -0.00404,-2.00002 0.00134,-2.79867c0.00533,-0.7907 0.01731,-1.94872 0.16533,-3.03466c0.09319,-0.68373 0.34203,-1.62286 0.48067,-2.205c0.24687,-1.03659 0.27422,-1.95678 0.51933,-2.96167c0.19685,-0.80698 0.32479,-1.61421 0.83334,-2.45833c0.4166,-0.6915 0.95807,-1.39709 1.71433,-1.54167c0.68231,-0.13043 1.68567,0 2.48566,0c1.00034,0 1.84167,-0.00081 2.57067,0c1.19866,0.00133 1.89266,-0.00184 2.61466,0.072c0.73479,0.07515 2.13816,0.47603 2.76801,0.78933c0.61053,0.30369 1.68277,0.9673 2.394,1.3905c0.68069,0.40502 1.62757,1.0719 2.13667,1.5175c0.55044,0.48179 1.09316,1.01347 1.632,1.548c1.002,0.99399 1.48185,1.42134 2.02233,2.13533c0.70093,0.92594 0.91077,1.50372 1.287,2.358c0.25989,0.59013 0.85565,1.79105 1.2,2.38934c0.34633,0.60174 1.02145,1.79387 1.28934,2.4c0.26265,0.59429 0.58412,1.49889 0.98534,2.70017c0.20166,0.60378 0.40953,1.21775 0.7,2.18517c0.31021,1.03316 0.5024,1.77121 0.7,2.5355c0.20245,0.78304 0.50093,1.97894 0.7,2.7805c0.20094,0.8091 0.40346,1.63386 0.60001,2.484c0.20357,0.88052 0.40079,1.7985 0.59999,2.74c0.1004,0.47458 0.29483,1.42313 0.60001,2.808c0.09818,0.44554 0.30316,1.30785 0.47917,2.14584c0.17113,0.8148 0.36578,2.01972 0.43683,2.8195c0.07016,0.78973 0.08537,1.56759 0.168,2.316c0.15518,1.40567 0.35588,2.04598 0.56334,2.98083c0.20006,0.9015 0.30716,1.80391 0.33183,2.4045c0.03696,0.89945 0.02083,2.1 0.02083,2.7c0,0.9 0,2.10001 0,2.7c0,0.9 0,2.1 0,2.7c0,0.9 0,1.8 0,2.70016c0,0.91051 0,1.857 0,2.8665c0,0.71867 0,1.46934 0,2.244c0,0.79066 0,1.98917 0,2.78484c0,0.78366 0.00339,1.91917 -0.00017,2.98299c-0.00221,0.66434 -0.01056,1.61341 -0.08516,2.22017c-0.14859,1.20855 -0.42987,2.0978 -0.56201,2.70133c-0.26138,1.19379 -0.32666,2.09966 -0.342,2.7c-0.03064,1.19965 -0.00085,2.09608 -0.03149,2.67916c-0.04398,0.83683 -0.28197,1.83104 -0.55383,2.74617c-0.27802,0.9358 -0.38495,1.73561 -0.41467,2.62133c-0.02147,0.63971 0.00256,1.67487 -0.02133,2.77867c-0.01624,0.74958 -0.03285,1.49686 -0.21066,2.553c-0.16537,0.98229 -0.40523,1.56387 -0.612,2.655c-0.17499,0.92336 -0.13547,1.56268 -0.168,2.56533c-0.03067,0.94532 -0.4382,1.80641 -0.71468,2.48534c-0.37103,0.91112 -0.25465,1.71789 -0.2885,2.7205c-0.02882,0.85381 -0.21069,1.81843 -0.42084,2.62083c-0.21571,0.82368 -0.48152,1.85023 -0.53867,2.53868c-0.06909,0.83252 0.06169,1.84583 -0.09315,2.57883c-0.22437,1.06207 -0.88196,1.75105 -0.94284,2.65717c-0.05363,0.7982 0,1.75717 0,2.65717c0,0.74733 0,1.7415 0,2.50684c0,1.07883 0,1.936 0,2.836c0,0.91467 0,1.7045 0,2.60001c0,0.96399 -0.23302,1.86946 0.0045,2.6955c0.22079,0.76787 0.82709,1.33496 1.01633,2.2045c0.18419,0.84628 0.76031,1.63635 1.06317,2.3c0.37695,0.82603 1.32007,1.40738 1.74933,2.2c0.36066,0.66595 1.02571,1.22804 1.25067,1.88667c0.32402,0.94868 0.65886,1.6738 1.31599,2.47733c0.55353,0.67683 1.22626,1.2536 1.70018,1.95734c0.446,0.66229 1.05251,1.52948 1.72516,2.05333c0.67706,0.5273 1.29469,1.10786 1.8,1.836c0.5457,0.78633 1.02184,1.39871 1.58533,1.98933c0.7758,0.81314 1.40738,1.38493 1.95334,2.036c0.36623,0.43674 0.955,1.44661 1.336,2.26401c0.25468,0.5464 0.80045,1.59473 1.336,2.27467c0.54749,0.69508 1.16611,1.39038 1.76851,2.17799c0.41779,0.54622 0.81462,1.13075 1.56216,1.88066c0.74968,0.75204 1.34789,1.14961 1.88066,1.61934c1.0614,0.93579 1.50481,1.50069 1.95267,2.06816c0.45761,0.57982 1.16191,1.48331 1.68401,2.08051c0.53645,0.61362 1.11046,1.21866 1.63199,1.884c0.79369,1.01256 1.2347,1.76805 1.46935,2.14c0.47368,0.75079 0.94336,1.52539 1.85733,2.5535c0.31194,0.35089 1.01515,0.99918 1.70667,1.6825c0.34898,0.34485 1.34315,1.45515 1.88133,2.308c0.27284,0.43237 0.74939,1.35104 0.98399,1.81601c0.46979,0.93109 0.93721,1.85687 1.47067,2.72932c0.26457,0.43271 0.82983,1.26859 1.41917,2.08084c0.29402,0.40523 1.1955,1.60449 1.4955,2.0045c0.59999,0.8 1.2,1.59999 1.5,2c0.89999,1.2 1.20657,1.59521 1.79549,2.4045c0.59048,0.81143 0.88702,1.21979 1.41917,2.08084c0.53761,0.8699 1.00407,1.79536 1.22816,2.26283c0.44444,0.92715 0.84747,1.84442 1.05583,2.2825c0.58665,1.23349 1.00035,1.95755 1.38051,2.6485c0.36705,0.66713 0.86964,1.68159 1.13683,2.40616c0.2747,0.74495 0.49602,1.52982 0.76801,2.31601c0.41718,1.20587 0.76497,2.00243 1.13683,2.79866c0.37544,0.80389 0.58122,1.19963 1.17917,2.40016c0.4005,0.80409 0.80182,1.6198 1,2.03584c0.40373,0.84758 1.00204,2.18462 1.2,2.64932c0.40421,0.9489 0.60172,1.4348 1,2.41917c0.20088,0.49648 0.79939,1.99557 1,2.49417c0.39879,0.99115 0.59639,1.48201 1,2.44417c0.19824,0.47258 0.79948,1.83582 1.19551,2.7c0.38333,0.83649 0.58488,1.23958 0.91916,2.05583c0.48796,1.19156 0.73479,1.99486 0.96451,2.7805c0.22276,0.76181 0.51358,1.85634 0.72083,2.53416c0.2899,0.94814 0.47369,1.53549 0.90001,2.56517c0.37698,0.91054 0.72513,1.5081 1.1955,2.5215c0.36719,0.79108 0.62486,1.59586 0.68367,2.39999c0.07295,0.99756 0.02083,1.79999 0.02083,2.60001c0,1.00449 0,1.85715 0,2.82132c0,0.80402 0.01147,1.62607 -0.00017,2.65735c-0.01041,0.92151 -0.04358,1.61453 -0.28517,2.60532c-0.2012,0.82513 -0.48869,1.64337 -0.91466,2.62933c-0.28214,0.65305 -0.76091,1.55276 -1.40001,2.276c-0.56689,0.64154 -1.13129,1.28952 -1.89999,1.88983c-0.65067,0.50812 -1.66667,0.52084 -2.52083,0.52084c-0.97917,0 -1.77917,0 -2.57918,0c-1,0 -1.79999,0 -2.80133,0c-0.61949,0 -1.76549,0.0119 -2.56033,-0.00449c-0.88536,-0.01825 -1.52902,-0.05316 -2.55949,-0.21701c-0.73038,-0.11615 -1.82176,-0.37433 -2.50951,-0.57983c-0.93271,-0.27869 -1.70653,-0.60217 -2.4175,-0.89868c-0.97134,-0.40512 -1.80862,-0.67825 -2.65634,-1.11432c-0.92039,-0.47345 -1.50699,-0.94287 -2.09549,-1.45767c-0.61229,-0.53561 -1.20226,-1.12842 -2.10001,-2.02817c-0.6015,-0.60284 -1.20235,-1.21844 -1.5,-1.53583c-0.6048,-0.64487 -1.5034,-1.68256 -2.09999,-2.41666c-0.30173,-0.37128 -0.89404,-1.12689 -1.5,-1.86334c-0.59419,-0.72217 -1.20041,-1.39597 -1.77867,-2.048c-0.55228,-0.62274 -1.10225,-1.20087 -1.6785,-2.13583c-0.53069,-0.86105 -0.6938,-1.49954 -1.30951,-2.53351c-0.44189,-0.7421 -0.8708,-1.13992 -1.34933,-1.85065c-0.65158,-0.96774 -0.97153,-1.8338 -1.17949,-2.42052c-0.42172,-1.18976 -0.59752,-1.80415 -0.94051,-2.73149c-0.24129,-0.65237 -0.67949,-1.68828 -1.02116,-2.42117c-0.35739,-0.76662 -0.74681,-1.54759 -0.94417,-1.94417c-0.59695,-1.19958 -1.00256,-1.99927 -1.37784,-2.81949c-0.18907,-0.41325 -0.57135,-1.25259 -0.9995,-2.60051c-0.29344,-0.92383 -0.39363,-1.40462 -0.64667,-2.35333c-0.12616,-0.47298 -0.52385,-1.86516 -0.85866,-2.73866c-0.32451,-0.84663 -0.51648,-1.24478 -0.86334,-2.02499c-0.48666,-1.0947 -0.90225,-2.10204 -1.13184,-2.74084c-0.32922,-0.91606 -0.62082,-1.82083 -0.92082,-2.72083c-0.2,-0.60001 -0.60001,-1.8 -0.8,-2.40001c-0.3,-0.89999 -0.5,-1.5 -0.90001,-2.7c-0.2,-0.59999 -0.5,-1.5 -0.7,-2.09999c-0.39999,-1.20001 -0.59999,-1.8 -0.89999,-2.70001c-0.3,-0.89999 -0.60001,-1.79999 -0.90001,-2.7c-0.2,-0.59999 -0.59999,-1.8 -0.79999,-2.39999c-0.3,-0.90001 -0.6091,-1.79706 -0.89551,-2.7c-0.18971,-0.59811 -0.44742,-1.4953 -0.58333,-2.40001c-0.13463,-0.89622 -0.16217,-1.81995 -0.4545,-2.7c-0.30623,-0.92189 -0.68469,-1.48085 -1.18401,-2.31467c-0.45844,-0.76555 -0.99529,-1.66748 -1.196,-2.36934c-0.28018,-0.97974 -0.21342,-1.80701 -0.37067,-2.54666c-0.24174,-1.13716 -0.75357,-1.63017 -1.31599,-2.60533c-0.32115,-0.55682 -0.82838,-1.54918 -1.164,-2.328c-0.3456,-0.80197 -0.3974,-1.58965 -0.44051,-2.49316c-0.04517,-0.94655 -0.29164,-1.737 -0.57016,-2.50684c-0.38976,-1.0773 -0.59396,-2.00192 -0.85068,-2.672c-0.29289,-0.76448 -0.45517,-1.77202 -0.62733,-2.50684c-0.21738,-0.92783 -0.78828,-1.82225 -0.922,-2.4825c-0.20366,-1.00571 -0.64764,-1.6315 -0.94133,-2.65866c-0.19223,-0.67233 -0.36972,-1.63176 -0.7885,-2.26868c-0.59154,-0.89964 -1.05798,-1.56313 -1.23834,-2.49016c-0.18313,-0.9413 -0.18095,-1.79158 -0.70449,-2.40001c-0.65479,-0.76094 -0.95459,-1.65254 -1.57417,-2.23582c-0.60439,-0.56898 -1.4785,-0.76868 -2.2785,-1.12134c-0.8,-0.35266 -1.66003,-0.73376 -2.5,-1.16667c-0.71348,-0.36772 -1.42533,-0.83333 -2.28533,-0.83333c-0.91334,0 -1.79867,0 -2.71467,0c-0.80133,0 -1.71741,0.02785 -2.47916,0.5c-0.76082,0.47157 -1.64685,0.69086 -2.22083,1.2045c-0.73534,0.65804 -1.34122,1.259 -1.9,1.9315c-0.54231,0.65266 -0.87551,1.41628 -1.56667,2.364c-0.45704,0.62669 -1.03326,1.21788 -1.34933,1.68533c-0.65283,0.96552 -0.97157,1.83244 -1.1795,2.41917c-0.42163,1.18973 -0.61046,1.7981 -0.9045,2.73149c-0.20415,0.64804 -0.50732,1.68336 -0.6955,2.41667c-0.19129,0.74541 -0.44621,1.85951 -0.58317,2.92599c-0.08562,0.66672 -0.13259,1.64807 -0.20533,2.30533c-0.11512,1.04021 -0.27048,1.77283 -0.416,2.53683c-0.14919,0.78333 -0.34899,1.97774 -0.41467,2.77783c-0.0648,0.78957 -0.08202,1.94417 -0.08533,3.03467c-0.00208,0.68268 -0.00389,1.66783 -0.05717,2.31934c-0.08312,1.01651 -0.22086,1.73549 -0.37883,2.494c-0.16696,0.80165 -0.3717,1.6386 -0.664,2.97484c-0.20446,0.93469 -0.30127,1.41423 -0.5,2.38916c-0.10064,0.4937 -0.39761,1.98982 -0.60133,2.98933c-0.10197,0.50024 -0.31502,1.50105 -0.43467,2c-0.24033,1.0022 -0.5224,2.00266 -0.848,3.00133c-0.1645,0.50455 -0.5192,1.51843 -0.69067,2.03467c-0.34783,1.0472 -0.67265,2.12724 -0.94,3.248c-0.13526,0.56702 -0.25524,1.1402 -0.47466,2.29066c-0.10912,0.57213 -0.20313,1.14249 -0.51083,2.80367c-0.09917,0.53539 -0.30499,1.57625 -0.42067,2.08c-0.22562,0.98253 -0.6386,2.37758 -0.80067,2.82051c-0.3214,0.87836 -0.69946,1.74193 -0.88917,2.1825c-0.39093,0.90788 -0.79453,1.87 -1.19067,2.908c-0.20587,0.53947 -0.61723,1.65236 -0.83466,2.21599c-0.21889,0.5674 -0.92563,2.25038 -1.18533,2.79317c-0.25769,0.53859 -0.80296,1.58607 -1.07867,2.10068c-0.27363,0.51073 -0.83021,1.51228 -1.33334,2.52083c-0.49542,0.99309 -0.71626,1.49707 -1.14583,2.5c-0.21354,0.49857 -0.6201,1.49976 -1.02084,2.5c-0.20013,0.49954 -0.59679,1.49049 -0.97916,2.45833c-0.18546,0.46942 -0.54048,1.36913 -0.9995,2.59868c-0.27161,0.72755 -0.50251,1.40916 -0.774,2.034c-0.5324,1.22531 -0.89002,1.81378 -1.26817,2.38815c-0.55944,0.84976 -1.17577,1.60211 -1.76849,2.28485c-0.55051,0.63409 -1.20214,1.42764 -1.68983,2.21515c-0.32832,0.53015 -0.84886,1.49384 -1.72083,1.93683c-0.84718,0.43039 -1.4792,0.94006 -2.34167,1.08401c-0.93211,0.15558 -1.705,0.00259 -2.625,0c-0.79501,-0.00226 -1.68157,-0.01871 -2.32,-0.08533c-1.03815,-0.10831 -2.15257,-0.34918 -2.924,-0.56134c-0.39558,-0.1088 -1.59615,-0.50345 -2.38667,-0.83865c-0.78847,-0.33432 -1.54371,-0.70279 -2.232,-1.116c-0.9466,-0.56827 -1.71153,-1.18475 -2.36166,-1.79868c-0.62158,-0.58698 -1.22329,-1.18274 -1.773,-1.80002c-0.6937,-0.77896 -1.24221,-1.57761 -1.47433,-2.19998c-0.36311,-0.97362 -0.35031,-1.84232 -0.52834,-2.76668c-0.15441,-0.80173 -0.42777,-1.67595 -0.54933,-2.31999c-0.19429,-1.02934 -0.24327,-1.76996 -0.27333,-2.924c-0.02059,-0.79047 -0.01067,-1.98933 -0.01067,-2.78934c0,-0.8 0,-1.60001 0,-2.8c0,-0.8 0,-1.60001 0,-2.39999c0,-0.8 0,-2 0,-2.8c0,-0.8 -0.001,-1.60001 0,-2.8c0.00067,-0.8 -0.0013,-1.60048 0.036,-2.39999c0.03736,-0.80066 0.17487,-2.00188 0.32117,-2.8c0.07351,-0.40102 0.24945,-1.20097 0.44417,-2c0.29275,-1.20132 0.49113,-1.99892 0.7195,-2.77917c0.22455,-0.76721 0.6271,-1.8632 0.96317,-2.53549c0.4891,-0.97842 0.89156,-1.56624 1.34933,-2.48534c0.44166,-0.88675 0.64465,-1.50972 1.05067,-2.39999c0.41523,-0.91049 0.76161,-1.50722 1.32667,-2.40001c0.3829,-0.60498 0.99206,-1.49837 1.38933,-2.10449c0.40277,-0.61452 1.00967,-1.5753 1.60017,-2.61684c0.20454,-0.36078 0.59982,-1.11017 1.03584,-1.85333c0.43867,-0.7477 1.17026,-1.81642 1.72117,-2.46817c0.54877,-0.64922 1.14597,-1.254 1.44417,-1.55583c0.89445,-0.90535 1.49941,-1.49626 2.07783,-2.10133c0.56887,-0.59508 1.36163,-1.49408 1.83684,-2.10001c0.46606,-0.59427 0.88177,-1.22214 1.85066,-2.39999c0.50329,-0.61183 1.05894,-1.19975 1.33334,-1.5c0.54787,-0.5995 1.58454,-1.79518 2.04283,-2.39549c0.44221,-0.57924 0.82706,-1.16776 1.65701,-2.18285c0.38405,-0.46973 0.96149,-1.07765 1.77866,-1.90016c0.49903,-0.50229 1.22308,-1.22258 1.77866,-1.77867c0.72583,-0.72649 1.3351,-1.16119 1.9095,-2.07616c0.4425,-0.70488 0.67577,-1.4805 1.33334,-2.16667c0.5871,-0.61263 1.30368,-1.29636 1.8,-1.8c0.70526,-0.71568 1.33694,-1.50188 1.67916,-2.2c0.39855,-0.81303 0.55103,-1.70786 0.94617,-2.39999c0.46472,-0.81404 1.06581,-1.41 1.54934,-2.03601c0.579,-0.74962 0.8663,-1.72908 1,-2.38933c0.20365,-1.00571 0.58419,-1.64914 1.12516,-2.49583c0.40893,-0.64003 0.91332,-1.51691 1.40017,-2.19966c0.51806,-0.72652 0.42247,-1.70399 0.5,-2.56317c0.08367,-0.92728 0.76117,-1.45367 0.97916,-2.416c0.19889,-0.87799 0.02084,-1.7955 0.02084,-2.66667l0,-0.89049l0,-0.89017l0,-0.886"></path>
                            <path id="chr_2" d="m344,115c0,-1.02084 0,-1.97884 0,-2.80983c0,-0.96084 0,-1.73183 0,-2.59283c0,-0.89734 0,-1.79733 0,-2.69733c0,-0.9 0,-1.8045 0,-2.736c0,-0.64934 0,-1.68517 0,-2.77867c0,-0.724 0,-1.416 0,-2.616c0,-0.7485 0,-1.936 0,-2.7265c0,-0.56367 0,-1.44417 0,-2.64283c0,-0.6 -0.0043,-1.50001 0.00015,-2.69984c0.0022,-0.59568 -0.04211,-1.4766 0.1665,-2.5335c0.18903,-0.95775 0.53775,-1.60746 0.71201,-2.58817c0.14349,-0.80749 0.13953,-1.68806 0.28799,-2.5785c0.15112,-0.90651 0.42801,-1.79565 0.54935,-2.4c0.23975,-1.19404 0.24438,-1.81196 0.33667,-2.75716c0.09872,-1.01111 0.21994,-1.72954 0.38333,-2.46817c0.24918,-1.12644 0.5618,-2.20864 0.76401,-2.892c0.19785,-0.66863 0.41901,-1.33897 0.79999,-2.816c0.10529,-0.40817 0.30203,-1.27605 0.5,-2.1875c0.20212,-0.93055 0.28897,-1.39898 0.60016,-2.7575c0.19736,-0.86163 0.39597,-1.68791 0.65701,-2.46c0.37442,-1.10748 0.7186,-1.78073 1.06818,-2.42567c0.34131,-0.62972 1.01913,-1.83088 1.29065,-2.436c0.40057,-0.89275 0.60806,-1.51371 1.05066,-2.4c0.45697,-0.91509 0.84552,-1.49311 1.18066,-2.1c0.65543,-1.18686 0.88763,-1.80101 1.27802,-2.664c0.25076,-0.55431 0.8356,-1.57109 1.39551,-2.21517c0.55698,-0.64071 1.38446,-1.41705 1.77914,-1.8315c0.60828,-0.63873 1.40024,-1.61063 1.8045,-2.142c0.41641,-0.54729 1.04916,-1.39525 1.81683,-2.12583c0.51117,-0.48646 1.08441,-0.89317 1.89468,-1.52017c0.76801,-0.59428 1.4549,-1.19415 2.0795,-1.80133c0.81342,-0.79075 1.21878,-1.18274 1.76849,-1.8c0.69373,-0.77896 1.16782,-1.62117 1.58865,-2.14283c0.66763,-0.82758 1.44641,-1.20994 2.22601,-1.75716c0.66904,-0.46963 1.51172,-1.08919 2.19598,-1.5c0.81192,-0.48746 1.474,-1.05248 2.15869,-1.56667c0.66901,-0.5024 1.47882,-1.02411 2.06934,-1.34933c0.88425,-0.48699 1.79407,-0.85689 2.39731,-1.07333c1.19382,-0.42835 1.80115,-0.60348 2.67917,-0.91067c0.83792,-0.29316 1.82669,-0.72057 2.5,-0.97917c1.0087,-0.3874 1.80957,-0.50171 2.53552,-0.5195c1.07431,-0.02633 1.896,-0.00133 2.80664,-0.00133c0.85333,0 1.81094,-0.09288 2.56268,0.08533c0.97745,0.23172 1.67645,0.65752 2.48001,1.31467c0.67682,0.55351 1.37711,1.13379 1.59998,2c0.22443,0.87231 0.03601,1.80133 0.03601,2.76667c0,0.83333 0,1.66667 0,2.586c0,0.9 0.00009,1.59467 0,2.622c-0.00009,0.94683 0.06509,1.66288 -0.16666,2.62533c-0.1976,0.82067 -0.63416,1.60026 -1.03336,2.4c-0.39948,0.80026 -0.91895,1.79113 -1.18933,2.4c-0.35049,0.7893 -0.53653,1.80743 -0.66333,2.65717c-0.14447,0.96812 -0.37112,1.7722 -0.55798,2.3175c-0.36566,1.06708 -0.67966,1.77277 -1.08936,2.62083c-0.39081,0.809 -0.79999,1.6045 -1.29999,2.6045c-0.29999,0.59999 -0.70001,1.4 -1.10001,2.2c-0.39999,0.8 -0.89999,1.8 -1.19998,2.4c-0.39999,0.80001 -0.80002,1.6 -1.20001,2.4c-0.29999,0.6 -0.85336,1.57729 -1.14282,2.4c-0.3436,0.97667 -0.26434,1.80429 -0.64117,2.716c-0.28003,0.67747 -0.56433,1.55829 -0.80002,2.4c-0.25024,0.89364 -0.92285,1.64654 -1.31598,2.30933c-0.40717,0.68647 -1.11389,1.49142 -1.43335,2.34133c-0.24069,0.64036 -0.18024,1.45667 -0.52383,2.43333c-0.28946,0.82272 -0.74283,1.6 -1.14282,2.4c-0.5,1 -0.79999,1.59999 -1.20001,2.4c-0.39999,0.79999 -0.80447,1.59776 -1.29999,2.6c-0.29822,0.60316 -0.70157,1.44467 -0.87866,2.42133c-0.14447,0.79672 -0.10889,1.63879 -0.34268,2.65717c-0.16003,0.69707 -0.54303,1.51021 -0.72598,2.5215c-0.14365,0.79402 -0.27939,1.59814 -0.55267,2.4c-0.34006,0.99774 -0.48199,1.78886 -0.5,2.83334c-0.01364,0.79189 0,1.65067 0,2.514c0,0.8955 0,1.85266 0,2.75266c0,0.9 0,1.7 0,2.6c0,0.94283 -0.00201,1.82084 0,2.68934c0.00226,0.96334 0.35809,1.66576 1.06134,2.18533c0.59979,0.44312 1.45599,0.42533 2.43866,0.42533c0.98666,0 1.75443,0.01016 2.66666,-0.16666c0.79932,-0.15494 1.65378,-0.45724 2.23468,-0.72c1.22238,-0.55293 1.81467,-0.95175 2.43466,-1.36c0.65509,-0.43137 1.68805,-1.17113 2.41666,-1.7105c0.37408,-0.27694 1.12875,-0.84513 1.86334,-1.44417c0.72409,-0.59048 1.40582,-1.19376 2.05865,-1.80933c0.65662,-0.61912 1.31,-1.27402 1.992,-1.95599c0.71936,-0.71932 1.09384,-1.09176 1.85419,-1.875c0.39154,-0.40334 0.79697,-0.80793 1.97916,-2.125c0.4046,-0.45075 0.80463,-0.91478 1.60001,-1.86934c0.40234,-0.48288 0.80112,-0.97375 1.19998,-1.4685c1.20343,-1.49274 1.60046,-1.99501 1.99869,-2.4955c0.39743,-0.49951 1.1824,-1.49854 1.56534,-2c0.38071,-0.49855 1.47217,-1.99805 1.81448,-2.49983c0.33899,-0.49692 0.98615,-1.4891 1.30066,-1.97933c0.31079,-0.48444 0.91153,-1.44078 1.52084,-2.35417c0.59229,-0.88787 0.89954,-1.31432 1.5,-2.16666c0.29974,-0.42549 0.90833,-1.2784 1.5,-2.16666c0.60727,-0.91167 0.90106,-1.38261 1.47916,-2.33334c0.28763,-0.473 0.84146,-1.41267 1.35419,-2.33334c0.48834,-0.87686 0.70523,-1.3099 1.16666,-2.125c0.44427,-0.78481 1.15091,-1.88727 1.68533,-2.55766c0.53516,-0.67132 0.83752,-0.98219 1.41916,-1.63667c0.89865,-1.01122 1.49658,-1.73154 1.80002,-2.1c0.61581,-0.74776 1.24435,-1.50223 1.58081,-1.86333c1.01187,-1.08598 1.39114,-1.41196 2.16135,-2.05867c0.78973,-0.66309 1.21594,-0.97678 2.51999,-1.992c0.45502,-0.35424 0.91776,-0.71682 1.85867,-1.45866c0.47488,-0.37441 0.94806,-0.75097 1.89066,-1.49067c0.91611,-0.71892 1.36124,-1.06416 2.22684,-1.72683c0.42056,-0.32197 1.24603,-0.94714 2.05585,-1.55583c0.80084,-0.60197 1.5997,-1.20041 1.99683,-1.49683c0.78305,-0.58448 1.54343,-1.15276 2.63785,-1.93783c0.34393,-0.24672 0.99661,-0.72529 1.95599,-1.32c0.60153,-0.37289 1.509,-0.8287 2.39999,-1.22134c0.57373,-0.25283 1.66367,-0.76447 2.35352,-1.278c0.6665,-0.49617 1.43866,-1.20938 2.06165,-1.59016c0.61459,-0.37566 1.81442,-0.86471 2.35718,-1.05717c0.88074,-0.31231 1.52274,-0.5079 2.55966,-0.8c0.72977,-0.20557 1.48056,-0.39307 2.59467,-0.7c1.06433,-0.29322 1.72241,-0.51163 2.67334,-0.764c0.60669,-0.16101 1.80966,-0.3648 2.41068,-0.4c0.89902,-0.05266 2.09998,-0.036 2.69998,-0.036c0.89999,0 2.10001,0 2.70001,0c0.89999,0 1.79999,0 2.69983,0c0.8895,0 1.74301,0 2.77884,0c0.92133,0 1.60532,0 2.604,0c0.79648,0 2,-0.13902 2.77567,0.02083c0.8877,0.18295 1.34265,0.7795 2.04166,1.47917c0.69968,0.70033 1.25803,1.22665 1.44284,2.15717c0.14597,0.73483 0.05716,1.74417 0.05716,2.57883c0,1.1855 0,1.72117 0,2.57467c0,1.18933 0,2.09383 0,2.71017c0,0.9645 0.0108,2.00105 -0.009,2.74083c-0.0209,0.78041 -0.05954,1.60851 -0.32434,2.905c-0.0936,0.45824 -0.36719,1.39088 -0.70831,2.33333c-0.17273,0.47724 -0.55252,1.41973 -0.95834,2.33333c-0.39313,0.88503 -0.8013,1.73557 -1.23599,2.556c-0.21652,0.40866 -0.92087,1.61711 -1.44934,2.40934c-0.26691,0.40015 -0.82745,1.18671 -1.41916,1.94417c-0.29291,0.37495 -1.19147,1.43946 -1.7955,2.1045c-0.59607,0.65627 -0.90802,0.97079 -1.5,1.63667c-0.91248,1.02637 -1.50049,1.76827 -1.80002,2.14133c-0.60098,0.74853 -1.19052,1.49869 -1.79999,2.208c-0.59222,0.68925 -1.20081,1.34068 -1.836,1.95599c-0.63956,0.61957 -1.32391,1.21885 -2.048,1.80933c-0.73459,0.59904 -1.11124,0.88687 -1.86334,1.44417c-1.07925,0.7997 -1.7695,1.2701 -2.41666,1.73183c-0.62747,0.4477 -1.24927,0.89185 -2.13599,1.64667c-0.61029,0.5195 -1.19412,1.08519 -1.80002,1.626c-0.89145,0.79569 -1.49011,1.27772 -2.08932,1.706c-0.84183,0.60169 -1.63687,1.09263 -2.35352,1.4895c-0.84451,0.46769 -1.67053,0.83282 -2.45715,1.278c-0.61099,0.34579 -1.396,0.96784 -2,1.47916c-0.9935,0.84107 -1.59125,1.25846 -2.20001,1.58933c-0.98682,0.53635 -1.61154,0.83433 -2.31433,1.11066c-1.11444,0.43819 -1.80106,0.60099 -2.52731,0.82084c-0.76846,0.23262 -1.96307,0.63415 -2.75702,0.96317c-0.79556,0.3297 -1.56253,0.69688 -2.31598,1.032c-1.05997,0.47145 -1.73044,0.72574 -2.3645,0.984c-0.92288,0.37589 -1.82434,0.81504 -2.42084,1.15266c-0.90543,0.5125 -1.48526,0.89319 -2.39984,1.32584c-0.87634,0.41456 -1.75406,0.6898 -2.77866,1.12167c-0.69348,0.2923 -1.34235,0.62251 -2.12152,1.22117c-0.62134,0.4774 -1.15982,1.08141 -2.19998,1.79333c-0.57935,0.39652 -1.19974,0.68453 -2.20001,1.18533c-0.79977,0.40042 -1.61121,0.77857 -2.37915,1.2c-0.86771,0.47616 -1.48605,1.0508 -2.34216,1.59983c-0.65448,0.41972 -1.53699,0.65353 -2.26801,0.9255c-0.86066,0.3202 -1.57834,0.97305 -2.24667,1.38534c-0.76648,0.47281 -1.564,1.07333 -2.164,1.58933c-0.60001,0.516 -1.38583,1.13938 -2.20001,1.58933c-0.78705,0.43496 -1.60007,0.81051 -2.19998,1.11066c-1.00146,0.50107 -1.62814,0.80785 -2.28534,1.18533c-1.01025,0.58027 -1.55249,0.93301 -2.362,1.51917c-0.79871,0.57835 -1.51962,1.16749 -2.10599,1.79549c-0.71655,0.76745 -1.08035,1.40884 -1.47202,2c-0.67285,1.01556 -1.19275,1.59627 -1.72198,2.14284c-0.75708,0.78189 -1.39569,1.30894 -1.94818,1.95717c-0.64209,0.75335 -1.06622,1.41427 -1.64734,2.05266c-0.4942,0.54291 -1.24265,0.77351 -1.92517,1.11534c-0.73062,0.36591 -1.63626,1.02621 -2.388,1.588c-0.41791,0.31232 -0.85397,0.64949 -1.776,1.37601c-0.48251,0.38019 -0.97266,0.77115 -2.46765,1.96782c-0.49789,0.39856 -0.99234,0.79443 -1.479,1.18951c-0.95129,0.77225 -1.41745,1.14134 -2.68799,2.24399c-0.39362,0.3416 -1.12271,1.01019 -1.46133,1.34134c-0.6579,0.64334 -1.28302,1.29759 -2.17218,2.347c-0.30591,0.36102 -0.91241,1.12349 -1.53583,1.925c-0.32217,0.41422 -1.33029,1.72028 -1.68552,2.175c-0.36346,0.46529 -1.11533,1.42433 -1.49933,1.92c-0.39069,0.50432 -0.78751,1.01712 -1.18048,1.545c-1.21054,1.62605 -1.60031,2.19753 -1.99869,2.77049c-0.40167,0.57768 -0.8009,1.16304 -1.19998,1.75351c-0.40091,0.5932 -1.59967,2.38956 -2.00018,2.98933c-0.40082,0.60022 -0.80295,1.20093 -1.21048,1.79999c-0.40881,0.60094 -0.82303,1.20164 -1.24652,1.8c-0.42584,0.60165 -1.30664,1.80215 -1.76398,2.40001c-0.9213,1.20438 -1.39508,1.79878 -1.858,2.39999c-0.46106,0.5988 -0.91849,1.19458 -1.3515,1.8c-0.8508,1.18954 -1.23175,1.79651 -1.59866,2.39867c-0.36169,0.5936 -1.03085,1.77951 -1.34933,2.36533c-0.31415,0.57787 -0.61176,1.15443 -1.22134,2.26933c-0.59329,1.08505 -0.89813,1.61058 -1.20001,2.13068c-0.59634,1.02745 -0.89917,1.53198 -1.79999,3.036c-0.29956,0.50017 -0.60278,0.99759 -1.18933,2c-0.29187,0.49881 -0.58716,0.99323 -1.12534,2c-0.53107,0.99347 -1.00418,1.99811 -1.22818,2.5c-0.44464,0.99632 -0.64993,1.4986 -1.05582,2.5c-0.40475,0.99861 -0.80133,2 -1.00134,2.5c-0.39999,1 -0.60001,1.5 -1,2.5c-0.39999,1 -0.60138,1.49962 -1,2.5045c-0.20068,0.50589 -0.60251,1.53055 -0.79999,2.05266c-0.6077,1.60663 -0.80209,2.16321 -0.99869,2.72551c-0.3952,1.13042 -0.58282,1.69574 -0.76532,2.256c-0.35522,1.09052 -0.67703,2.12579 -0.95065,3.092c-0.25659,0.90599 -0.37192,1.34669 -0.58084,2.2175c-0.31573,1.31587 -0.50549,2.23695 -0.60449,2.709c-0.20203,0.96326 -0.40027,1.9473 -0.5,2.44417c-0.30072,1.4982 -0.39679,1.99667 -0.60001,2.97734c-0.19687,0.94998 -0.28339,1.41101 -0.60001,2.688c-0.19006,0.76653 -0.40121,1.47163 -0.70001,2.48149c-0.30124,1.01814 -0.50467,1.73729 -0.69998,2.4985c-0.20493,0.79875 -0.4043,1.63904 -0.70016,2.97484c-0.10266,0.46349 -0.30396,1.41603 -0.53583,2.38916c-0.11786,0.49464 -0.52597,1.99081 -0.68533,2.4895c-0.3219,1.0074 -0.50833,1.50903 -0.85333,2.53583c-0.17487,0.52043 -0.66336,2.12947 -0.80402,2.68533c-0.14206,0.56149 -0.396,1.70401 -0.52133,2.27867c-0.12534,0.57466 -0.36301,1.72153 -0.66666,2.83333c-0.29691,1.08713 -0.48383,1.61024 -0.65869,2.13068c-0.345,1.02679 -0.55057,1.52255 -1.008,3.036c-0.30017,0.99316 -0.41556,1.49406 -0.64581,2.47917c-0.11316,0.48415 -0.31201,1.43718 -0.62085,2.79916c-0.19525,0.86108 -0.29507,1.27625 -0.5,2.075c-0.19531,0.76122 -0.59586,2.16263 -0.79999,2.81067c-0.29404,0.93338 -0.59421,1.83798 -0.90451,2.74051c-0.21237,0.61771 -0.43259,1.26114 -1.01682,2.61682c-0.15732,0.36502 -0.5058,1.10318 -0.85333,1.85333c-0.34308,0.74054 -0.94394,2.13715 -1.18936,2.77866c-0.229,0.5986 -0.60693,1.69321 -0.94049,2.38953c-0.41705,0.87064 -1.01166,1.66153 -1.52084,2.25714c-0.51678,0.60449 -1.16232,1.39401 -1.65381,2.17917c-0.53082,0.84802 -1.00046,1.52036 -1.62085,2.12067c-0.58154,0.56271 -1.49994,0.89484 -2.42084,0.90018c-0.80048,0.00464 -1.65524,0.05057 -2.59515,-0.00134c-0.77142,-0.0426 -1.73465,-0.39038 -2.39999,-0.79999c-0.89661,-0.552 -1.38843,-1.0643 -2.13669,-1.86032c-0.53699,-0.57126 -1.20105,-1.49445 -1.46198,-2.02234c-0.39603,-0.80118 -0.66318,-1.64066 -1.08536,-2.63068c-0.28598,-0.67062 -0.82068,-1.67566 -1.09549,-2.28983c-0.37604,-0.84039 -0.5929,-1.81573 -0.67917,-2.62085c-0.08795,-0.82098 -0.23578,-1.85663 -0.45065,-2.54932c-0.21136,-0.6814 -0.5043,-1.91089 -0.53867,-2.47202c-0.05594,-0.91382 -0.03601,-1.91998 -0.03601,-2.63867c0,-1.13548 0,-1.91916 0,-3.11465c0,-0.8 0,-1.60001 0,-2.39999c0,-0.8 0,-2 0,-2.8c0,-0.8 0,-1.60001 0,-2.8c0,-0.79866 0,-1.58934 0,-2.364c0,-0.75067 0,-1.8145 0,-2.8c0,-0.62534 0,-1.836 0,-2.43733c0,-0.91949 -0.00299,-1.884 0.00015,-2.92c0.00223,-0.73135 -0.00381,-1.48027 0.05701,-2.226c0.11713,-1.43579 0.29449,-2.09573 0.44284,-2.73184c0.21454,-0.91985 0.41931,-2.12051 0.47467,-2.72083c0.08298,-0.89957 0.16245,-1.80676 0.38251,-2.7c0.22369,-0.90787 0.44357,-1.5002 0.84283,-2.70001c0.19968,-0.60008 0.5,-1.5 0.79999,-2.39999c0.29999,-0.89999 0.5,-1.5 0.79999,-2.39999c0.30002,-0.90001 0.60001,-1.8 0.80002,-2.40001c0.39999,-1.2 0.69574,-2.09697 0.89999,-2.67917c0.29388,-0.83766 0.68961,-1.84029 1,-2.5c0.48392,-1.02855 0.78949,-1.63348 1.25717,-2.478c0.39688,-0.71664 0.89594,-1.49522 1.24283,-2.04283c0.51663,-0.81558 1.12463,-1.80928 1.47916,-2.47917c0.53607,-1.01292 0.82083,-1.62082 1.22086,-2.42082c0.39999,-0.79999 0.80051,-1.59975 1.09998,-2.20001c0.50089,-1.00406 0.81699,-1.62888 1.20001,-2.56667c0.31873,-0.78033 0.57437,-1.6412 1.16666,-2.66666c0.28836,-0.49925 0.85333,-1.12973 1.34933,-1.76534c0.76285,-0.97755 1.07547,-1.62532 1.65067,-2.40134c0.45981,-0.62033 1.06207,-1.15883 1.81183,-2.2c0.41675,-0.57875 0.71014,-1.20351 1.22153,-2.19867c0.39471,-0.76808 0.87433,-1.60176 1.39999,-2.19067c0.66431,-0.74423 1.15323,-1.28616 1.37466,-2.23599c0.23297,-0.99931 0.58563,-1.64822 1.12534,-2.49617c0.43259,-0.67964 0.70023,-1.27971 1.19998,-2.2785c0.30017,-0.59991 0.79205,-1.60435 1.11069,-2.21066c0.46295,-0.88092 0.905,-1.61455 1.42532,-2.42532c0.37051,-0.57735 0.75925,-1.16943 1.56534,-2.36401c0.40628,-0.60207 0.82199,-1.2097 1.284,-1.79999c0.71667,-0.91563 1.27222,-1.50018 1.81467,-2.12083c0.55234,-0.63199 1.34662,-1.63712 1.81467,-2.36317c0.23492,-0.36441 0.66464,-1.11592 1.08084,-1.86333c0.60471,-1.08589 0.99796,-1.77139 1.40448,-2.41667c0.59061,-0.9375 1.20096,-1.83537 1.59869,-2.43466c0.3916,-0.59004 0.95682,-1.43645 1.43466,-2.23466c0.42075,-0.70288 0.86731,-1.55959 1.26648,-2.3665c0.38431,-0.7769 0.8121,-1.62238 0.89569,-2.49567c0.08621,-0.90041 0.31033,-1.86267 -0.0062,-2.69383c-0.28799,-0.75626 -1.30392,-0.84995 -2.01016,-1.11067c-0.97034,-0.35822 -1.78934,-0.48076 -2.66449,-0.49867c-1.00244,-0.02052 -1.87717,-0.00133 -2.53351,-0.00133c-1.1145,0 -1.95187,-0.00375 -2.85315,0.01067c-0.46973,0.00751 -1.43924,0.01953 -2.92801,0.15601c-0.50107,0.04593 -1.50565,0.18835 -2.01068,0.28c-1.02786,0.18652 -1.55188,0.29414 -3.15598,0.72c-0.55664,0.14778 -1.68738,0.4908 -2.25867,0.67999c-0.57608,0.19078 -1.15787,0.37851 -2.29068,0.84c-1.12622,0.4588 -1.6676,0.72459 -2.72684,1.26601c-0.52371,0.2677 -1.54599,0.82292 -2.05267,1.09467c-0.50266,0.26958 -2.00226,1.03592 -2.50449,1.26733c-0.9957,0.4588 -1.49823,0.66878 -2.5,1.08083c-0.49911,0.20531 -2.00046,0.80438 -2.50266,1.00584c-1.01892,0.40874 -1.53931,0.61874 -2.06934,0.83466c-1.0993,0.4478 -2.26169,0.92989 -2.87085,1.1855c-0.62683,0.26303 -1.26849,0.53223 -2.578,1.09933c-0.66525,0.2881 -1.33209,0.58122 -1.99649,0.88051c-1.31821,0.5938 -1.96188,0.89751 -2.59866,1.19867c-0.6319,0.29884 -1.88382,0.90038 -2.50934,1.2c-0.62714,0.30038 -1.25626,0.60516 -2.54132,1.2c-0.65344,0.30247 -1.31647,0.60065 -1.98602,0.89549c-0.67178,0.29584 -2.02103,0.87082 -2.69464,1.14734c-0.66772,0.27411 -1.98499,0.79327 -2.63101,1.03551c-0.634,0.23773 -1.87399,0.67801 -2.48001,0.87999c-0.59402,0.19798 -1.17422,0.40021 -2.875,0.875c-0.54321,0.15164 -1.6044,0.41193 -2.12,0.53067c-0.99893,0.23006 -1.48608,0.32652 -2.88,0.636c-0.88071,0.19553 -1.71835,0.39703 -2.51999,0.564c-0.75851,0.15797 -1.82382,0.33029 -2.4895,0.37885c-0.94513,0.06894 -2.15717,0.05716 -2.75583,0.05716c-0.88049,0 -1.716,0 -2.716,0c-0.87466,0 -1.888,0.04736 -2.68533,-0.02083c-1.0069,-0.08612 -1.81482,-0.36847 -2.51468,-0.7805c-0.76373,-0.44963 -1.37256,-0.99371 -2.31065,-1.60933c-0.52382,-0.34375 -1.3019,-1.06487 -1.76401,-1.76401c-0.52512,-0.79446 -1.05434,-1.37442 -1.70616,-2.12534c-0.57866,-0.66664 -0.85167,-1.466 -1.27916,-2.336c-0.40312,-0.82042 -0.33794,-1.57718 -0.70934,-2.48c-0.33524,-0.81491 -0.69682,-1.71942 -0.71468,-2.48534c-0.02327,-0.9984 -0.00133,-1.79866 -0.00133,-2.79866c0,-0.8 -0.00006,-1.79984 0,-2.57918c0.00009,-0.9995 -0.00589,-1.79703 0.35716,-2.67349c0.29723,-0.71754 0.82362,-1.48187 1.40001,-2.10001c0.70511,-0.75616 1.27039,-1.28423 1.92949,-1.84866c0.72002,-0.61659 1.26401,-1.01523 2.16,-1.59866c0.63246,-0.41183 1.66939,-1.01013 2.41051,-1.39551c0.75824,-0.39429 1.14529,-0.58762 1.94415,-0.91916c1.18716,-0.49271 1.99869,-0.73267 2.79868,-0.98534c0.79999,-0.25267 2.00421,-0.65633 2.79733,-0.98399c0.39584,-0.16353 1.1637,-0.51045 1.88834,-0.86334c1.27612,-0.62143 2.04233,-1.01666 2.70532,-1.34816c0.80902,-0.40451 1.60751,-0.80751 2.60901,-1.3045c0.59912,-0.2973 1.38428,-0.7088 2.39999,-0.916c0.79117,-0.16139 1.80215,-0.29912 2.39999,-0.50934c0.80298,-0.28236 1.81033,-0.44551 2.54285,-0.62733c1.06476,-0.26431 1.81537,-0.84413 2.67801,-0.9265c0.97495,-0.09308 1.79285,-0.00026 2.43631,-0.07799c0.97116,-0.11731 1.77402,-0.38799 2.64734,-0.59018c0.89105,-0.2063 1.80029,-0.29694 2.73151,-0.37799c0.64951,-0.05654 1.68668,-0.18798 2.41217,-0.32734c1.09286,-0.20993 1.77783,-0.39372 2.98117,-0.56333c0.72665,-0.10242 1.77621,-0.23406 2.57065,-0.50933c0.60449,-0.20946 1.59909,-0.52184 2.436,-0.564c0.92953,-0.04684 1.71667,-0.01067 2.53867,-0.01067c1.05869,0 2.00314,0.02295 2.67801,-0.05717c0.9758,-0.11584 1.7785,-0.38796 2.65182,-0.59016c0.89105,-0.20631 1.79285,-0.28185 2.65952,-0.378c1.07205,-0.11893 1.76672,-0.34882 2.57434,-0.622c1.01685,-0.34396 1.862,-0.29561 2.82831,-0.35267c0.72379,-0.04274 1.65338,-0.47143 2.43466,-0.716c0.78543,-0.24587 1.59882,-0.26583 2.59869,-0.28416c0.80011,-0.01467 1.815,-0.12543 2.59998,-0.43583c0.81787,-0.32341 1.57925,-0.81804 2.436,-1.12801c0.91107,-0.32961 1.71567,-0.35368 2.50269,-0.46133c1.12488,-0.15385 1.77872,-0.72145 2.66132,-0.93867c0.8746,-0.21526 1.80008,-0.02534 2.69983,-0.036c0.87952,-0.01042 0.89804,-0.95534 1.66183,-1c0.88953,-0.05202 1.80988,0.02803 2.59167,-0.42533c0.80344,-0.46592 1.64331,-0.57354 2.22134,-1.14934c0.64807,-0.64556 1.28793,-1.17049 1.90448,-1.9045c0.58319,-0.69426 0.84036,-0.87666 1.5,-1.52084c0.62869,-0.61395 0.75891,-1.38367 1.02084,-2.2215c0.26208,-0.83832 0,-1.78917 0,-2.6945l0,-0.8l0,-0.93133l0.35266,-0.70534"></path>
                            </g>
                       </svg>
                </div>
            </ul>
        </div>
                  <div class="chromoMap_toolarea">
                   <div id="chromoMap_panel001" data-tab-content class="visible">
                    <div class="panel">

                    <div class="area">
                   
                    <h2 style="padding-right:5px;"> chromoMap JS (Alpha) </h2>
                    <br>
                    <div class="divider"></div>
                    
                    <a href="https://lakshay-anand.github.io/chromoMap/index.html" style="padding-left:5px;color:black;text-decoration:none;" target="_blank" rel="noopener noreferrer"> Copyright <br> 2024 <br> Lakshay Anand. </a>
                    </div>
               
                        <div class="divider"></div>
                        <div class="area">
                            <p>Zoom </p>
                            <button id="zm">plus</button>
                            <button id="zmminus">minus</button>
                        </div>
                       
                        <div class="divider"></div>
    
                        <div class="area">
                            <center><p> Genomic <br> Range </p></center>
                            <input type="text" class="text">
                            <button>Go</button>
                        </div>
    
                        <div class="divider"></div>
    
                        <div class="area">
                            <p>Gene or <br> feature </p> <br>
                            <input type="text" class="text">
                            <button>Search</button>
                        </div>
                    </div>
                   </div>
                   <div id="chromoMap_panel002" data-tab-content>
                    <div class="panel">
                        
                    </div>
                   </div>
                   <div id="chromoMap_panel003" data-tab-content>
                    <div class="panel">
                        <div class="area">
                            <p>add marker</p>
                            <button id="add_marker">+</button>
                        </div>
                        <div class="divider"></div>
                    </div>
                   </div>
                   <div id="chromoMap_panel004" data-tab-content>
                    <div class="panel">
                    <div class="area">
                     <h3 style="padding-right:5px;"> Export Plot </h3>
                    
                    
                    <input type="text" id="flnme" value="my_chromoMap_Plot">
                    <button id="save_png">PNG</button>
                    <button id="save_svg">SVG</button>
                    
                    </div>
                    
                    <div class="area">
                    <p> Export <br> Settings </p>
                    <label> background: <label>
                    <input type="color" id="bgcoll" value="#ffffff"></input>
                    <button id="cmap_bg_col">change</button><br>
                    <label> Scale(1-6) </label>
                    <input type="range" min="1" max="6" value="3" id="imgscl" />

                    
                    </div>
                    <div class="divider"></div>
                    <div class="area">
                     <h3 style="padding-right:5px;"> Export Data </h3>
                     <button>CSV</button>
                     <button>JSON</button>
                     </div>

                    </div>
                   </div>
                  </div>
    
        </div>



`;
var style = document.createElement('style');
style.appendChild(document.createTextNode(css));

document.getElementsByTagName('head')[0].appendChild(style);


document.getElementById("chromoMap_Area").innerHTML += ch_tool_bar;

const list = document.querySelectorAll(".list");
function activeLink(){
    list.forEach(
        (item) => 
        item.classList.remove("active")   
    );
    this.classList.add('active')
}
list.forEach( (item) => 
item.addEventListener('click',activeLink))

const tabs = document.querySelectorAll('[data-tab-target]');
const tabContents = document.querySelectorAll('[data-tab-content]');
tabs.forEach(tab => {
tab.addEventListener('click', () => {
    const target = document.querySelector(tab.dataset.tabTarget);
    tabContents.forEach(tabcontent => tabcontent.classList.remove('visible'))
    target.classList.add('visible');
})
})

document.getElementById('cmap_canvas').innerHTML = '';



d3.select("#cmap_canvas").append("div")
          .attr("id", div_id);
        

d3.select('body').style("overflow-y","scroll")
.style("overflow-x","scroll");

var chmap_div = d3.select("#"+div_id);
chmap_div.style("background-color","white");



var svg = chmap_div.append("svg")
  .attr("width",width)
  .attr("height",height)
  .attr("id",div_id+"-svg"); 


const mrkr = document.getElementById("add_marker");

var deltaX, deltaY;
var dragHandler = d3.drag()
                  .on("start", 
                  function() {
                    d3.select(this).style("fill", "black");
                    var current = d3.select(this);
                    deltaX = current.attr("cx") - d3.event.x;
                    deltaY = current.attr("cy") - d3.event.y;
                  }
                  )
                  .on("drag", 
                  function() {
                    d3.select(this)
                    .attr("cx", d3.event.x + deltaX)
                    .attr("cy", d3.event.y + deltaY);                    }
                  )
                  .on("end", 
                  function() {
                    d3.select(this).style("fill", "red");
                  }
                  );

                  mrkr.onclick = function(){


                    svg.append("circle")
                      .style("stroke", "gray")
                      .style("fill", "red")
                      .attr("r", 15)
                      .attr("cx", 50)
                      .attr("cy", 20)
                      .attr("id", "marker_1")
                      .call(dragHandler); }






let btn = document.getElementById("save_png");

btn.onclick = function(){
saveSvgAsPng(document.getElementById(div_id+"-svg"), document.getElementById("flnme").value + ".png", {scale: document.getElementById("imgscl").value,backgroundColor:document.getElementById("bgcoll").value,encoderOptions:1});
}


let btn2 = document.getElementById("save_svg");

btn2.onclick = function(){
saveSvg2(div_id+'-svg', document.getElementById("flnme").value + ".svg");
}

let btnbg = document.getElementById("cmap_bg_col");

btnbg.onclick = function(){

  chmap_div.style("background-color",document.getElementById("bgcoll").value);
  d3.selectAll(".overlines2").attr("fill",document.getElementById("bgcoll").value);

}

//btn.style = "background-color: #4CAF50; color: white;padding: 15px 32px;text-align: center;text-decoration: none;display: inline-block;font-size: 16px;width:100%;border:solid;border-radius:10px;"
//cmap_div_ele = document.getElementById(div_id);
//cmap_div_ele.appendChild(dwdiv);

///

var deltaX, deltaY;
var dragHandler = d3.drag()
                    .on("start", 
                    function() {
                      d3.select(this).style("fill", "black");
                      var current = d3.select(this);
                      deltaX = current.attr("cx") - d3.event.x;
                      deltaY = current.attr("cy") - d3.event.y;
                    }
                    )
                    .on("drag", 
                    function() {
                      d3.select(this)
                      .attr("cx", d3.event.x + deltaX)
                      .attr("cy", d3.event.y + deltaY);                    }
                    )
                    .on("end", 
                    function() {
                      d3.select(this).style("fill", "red");
                    }
                    );
 
                    
 // UI component as JS class 
 class UI {

  constructor(){
      this.ui_obj = undefined;
  }
  createButton(name){
      this.ui_obj = document.createElement("button");
      this.ui_obj.innerHTML = name;
      return this;
  }

  createDiv(){
    this.ui_obj = document.createElement("div");
    return this;
  }

  set(property,value){
     this.ui_obj[property] = value;
     return this;
  }

  get(){
    return this.ui_obj;
  }

 } 

 // UI end 
 
 var ui = new UI();

 let btn3 = ui.createButton("mybutton")
         .set("id","my_btn1")
         .set("class","btnnn")
         .get();



}

if(!v_align){
/* Adding title */
svg.append("text")
.attr("x",width/2)
.attr("y",(top_margin/2))
.attr("font-family", "sans-serif")
.attr("font-size", title_font_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(ttl);
} else {
svg.append("text")
.attr("x",width/2)
.attr("y",(left_margin/3))
.attr("font-family", "sans-serif")
.attr("font-size", "12px")
.attr("fill", "black")
.style("text-anchor", "middle")
.text(ttl);

}


if(!v_align){
/*adding axes  */
x_scale_pos=( top_margin+(chr_spacing)*nLoci[0].length+10)- y_chr_scale;


nnmax= win_scale;

if(cnt){
x_scale_width= (nnmax)*loci_width ;
} else {
x_scale_width= (nnmax)*loci_width ;
}


// Create scale
cc_ary = d3.range(0,x_scale_width,loci_width)
cc_ary.push(x_scale_width)
console.log(cc_ary)
var scale = d3.scaleLinear()
 .domain(ch_domain)
 .range(cc_ary);
//var scale_suffix = "bp"
//console.log(scale(700));
//console.log((parseInt(ch_domain[1])- parseInt(ch_domain[0])));
var show_pos_scale = false;
if(show_pos_scale){
var p_scale = d3.scaleLinear()
   .domain([1,((parseInt(ch_domain[1])- parseInt(ch_domain[0])))]).nice()
   .range([0, x_scale_width]);

   var x_axis2 = d3.axisBottom()
                .scale(p_scale)
                .ticks(5)
                .tickFormat(function(d){
                  if(d>=1000000000){ return (d/1000000000)+"Gb"
                  } else if(d>=1000000){

                  return (d/1000000)+"Mb";
                } else if (d>=1000) {
                  return (d/1000)+"kb";
                } else { return d+"bp";}});
}

//var tickValues= [].concat(scale.domain()[0], scale.ticks(5), scale.domain()[1]);
// Add scales to axis
var x_axis = d3.axisBottom()
.scale(scale)
.ticks(scale_ticks)                
.tickFormat(function(d){
  if(d>=1000000000){ return (d/1000000000)+"G"+scale_suffix
  } else if(d>=1000000){

  return (d/1000000)+"M"+scale_suffix;
} else if (d>=1000) {
  return (d/1000)+"k"+scale_suffix;
} else { return d+scale_suffix;}});

//Append group and insert axis
svg.append("g")
.attr("transform", "translate("+(left_margin)+"," + x_scale_pos + ")")
.attr("id","main_axis")
.call(x_axis);


svg.select("#main_axis").style("font-size","12px")
if(show_pos_scale){
svg.append("g")
.attr("transform", "translate("+(left_margin)+"," + (x_scale_pos+30) + ")")
.attr("id","main_axis2")
.call(x_axis2);

svg.select("#main_axis2").style("font-size","15px")
}

/* axis title */
/*
svg.append("text")
.attr("x",x_scale_width/2)
.attr("y",x_scale_pos + 45 + 15)
.attr("font-family", "sans-serif")
.attr("font-size", "15px")
.attr("fill", "black")
.style("text-anchor", "middle")
.text("Length");
*/

var grid_h = 0 - (x_scale_pos - top_margin + 15);
//var grid_array = [0,20,56,80,90];
/*  adding vertical grid lines*/
//console.log(scale(150))
//var grid_text = ["AUREA",""];
//var liness = [];
if(typeof grid_array == 'number'){
grid_array = [grid_array];
grid_text = [grid_text];
}
//console.log(grid_array.length)

if(vertical_grid){
for(i=0;i<grid_array.length;i++){

svg.select("#main_axis").append("line")
.attr("x1", scale(grid_array[i]))
.attr("x2", scale(grid_array[i]))
.attr("y1",  0)
.attr("y2",  grid_h)
.attr("stroke", grid_color)
.attr("stroke-dasharray","5,5")
.attr("stroke-width",2)
.attr("id","grid_id_"+i)
.attr("class","overlines");

let line_x = parseFloat(d3.select("#grid_id_"+i).attr("x1"));

if(grid_text != null){
svg.append("text")
.attr("x",line_x+left_margin+4)
.attr("y",grid_text_y)
.attr('text-anchor', 'start')
.attr("font-size",grid_text_size)
.attr("font-family", "sans-serif")
.attr("class", "myLabel")
.text(grid_text[i]);
}

}}
/* */
// window guides
if(guides){


  var geno_tip  = d3.select("body").append("div")
.attr("class", div_id+"-"+"genotooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 180px;height: 140px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;

  var chDataReducedPerRegion = d3.nest()
  .key(function(d) { return d.loci_start; })
  .rollup(function(v) { var label =''; var tot = '<font color="grey"> Count: </font> '+v.length+"  </p></div><hr><div  style='float:top;height:80%;overflow-wrap: break-word;' > <font size='1'>";for(var i=0; i< v.length;i++){ label = label +"<a href="+v[i].hlink+" style='text-decoration:none;cursor:pointer;color:grey;pointer-events: all'>"+ v[i].name+"</a>" +" ,";    }; return tot+label.slice(0,-1);})
  .entries(chData.flat(1));




var chDataReducedPerRegionALL = [];

for(i=0;i<ch_domain.length;i++){
   
  var t_arr = chDataReducedPerRegion.filter(object => object.key == String(ch_domain[i]));
  if(t_arr.length != 0){
    chDataReducedPerRegionALL[i] = {
      region_start: ch_domain[i] ,
      region_end: ch_domain[i+1]-1,
      value: chDataReducedPerRegion.filter(object => object.key == String(ch_domain[i]))[0].value}
  } else {
  chDataReducedPerRegionALL[i] = {
    region_start: ch_domain[i] ,
    region_end: ch_domain[i+1]-1,
    value: 'NA' }
 
  }
  }

print(chDataReducedPerRegionALL)

for(i=0;i<ch_domain.length;i++){



svg.select("#main_axis").append('rect')
.datum(chDataReducedPerRegionALL[i])
.attr("x", scale(ch_domain[i]))
.attr("height",  0 - grid_h)
.attr("y",  grid_h)
.attr("fill", "white")
.attr("width",chr_length)
.attr("id","guide_id_"+i)
.attr("class","overlines2");

}

d3.selectAll(".overlines2").on("mouseover",function(d){
  d3.select(this).attr("fill",guides_color);
  geno_tip.transition().duration(200).style("opacity", .9)
               .style("background-color","#F5F5F5")
               .style("width",180+'px');

               geno_tip.html("<div  style='float:top;height:10%;color:black;font-size:12px;'><p>"+d.region_start+" - "+d.region_end+"</p></div><hr><div  style='float:top;height:10%;'><p ><b>"+d.value+"</font></div>")
               .style("left", (d3.event.pageX) + "px")
               .style("top", (d3.event.pageY - 18) + "px");
})
.on("mouseout",function(d){
  d3.select(this).attr("fill",chmap_div.style("background-color"));
  geno_tip.transition() .delay(1000).duration(500)	.style("opacity", 0);
});


}

/* end of guides */

} else {

/*adding axes  */
x_scale_pos=( top_margin+(chr_spacing)*nLoci[0].length+10);



nnmax= 100;

if(cnt){
x_scale_width= (nnmax)*loci_width + 3*arc_radius;
} else {
x_scale_width= (nnmax)*loci_width + arc_radius;
}
// Create scale
var scale = d3.scaleLinear()
.domain(ch_domain).nice()
.range([0, x_scale_width]);

// Add scales to axis
var x_axis = d3.axisRight()
.scale(scale)
.ticks(5)
.tickFormat(function(d){
if(d>=1000000000){ return (d/1000000000)+"Gb"
} else if(d>=1000000){

return (d/1000000)+"Mb";
} else if (d>=1000) {
return (d/1000)+"kb";
} else { return d+"bp";}});

//Append group and insert axis
svg.append("g")
.attr("transform", "translate("+x_scale_pos+"," + (left_margin) + ")")
.call(x_axis).selectAll("text").attr("transform", "rotate(90)").attr("y", -12).style("text-anchor", "middle");


/* axis title */
svg.append("text")
.attr("transform", "translate("+(x_scale_pos + 35)+"," + (x_scale_width/2) + ") rotate(90)")
.attr("font-family", "sans-serif")
.attr("font-size", "12px")
.attr("fill", "black")
.style("text-anchor", "middle")
.text(" Length (bp)");




}


// scatter color function and scatter legend

//console.log(uniq_cates)
scatter_color_map_fn = []
/* ################################
scatter_color_map_fn = d3.scaleOrdinal()
.domain(uniq_cates)
.range(scatter_col);  */


if(scatter_mapping){


w2=uniq_cates.length*20;
// Create scale

var scalescatter = d3.scaleBand()
.domain(uniq_cates)
.range([2,w2]);

// Add scales to axis
var legendscatter = d3.axisRight()
.scale(scalescatter).tickSizeOuter(0);


//Append group and insert axis
svg.append("g")
.attr("transform", "translate("+((width-35)-scatter_lg_x)+"," + ((height-100)-scatter_lg_y )+ ")")
.attr("class","scatter_legend")
.call(legendscatter);

svg.select(".scatter_legend").style("font-size","15px");

//ading the label for the scatter
svg.append("text")
 .attr("x",(width-35)-scatter_lg_x)
 .attr("y",((height-100)-scatter_lg_y - 15 ) )
 .attr("font-family", "sans-serif")
 .attr("font-size", 20)
 .attr("fill", "black")
 .style("text-anchor", "start")
 .text(cat_legend_lab)
 .attr("class","sctrlbl");

// svg.append("g")
// .attr("transform", "translate("+width/2+"," +height/2+ ")")
// .attr("id","scatter_legend")
// .call(legendscatter);

rec_h2=w2/uniq_cates.length;
 svg.selectAll(".rects")
 .data(scatter_col)
 .enter()
 .append("rect")
 .attr("y", function(d,i){return (((height-99)-scatter_lg_y)+i*rec_h2);})
 .attr("height", rec_h2)
 .attr("x", function(d,i){return ((width-45)-scatter_lg_x);})
 .attr("width", 10)
 .attr("fill", function(d){return d;});


}

/*  Rendering chromoMap  */
u=1;
while(u<=parseInt(ploidy_n)){

renderChromoMap(chData[u-1],chr_spacing,ploidy=u,ch_width,svg,nLoci[u-1],cnt,ch_gap,chr_col[u-1],
       heatmap,left_margin,lg_x,lg_y,heat_scale,labels,width,height,rng[u-1],heat_col[u-1],
       an_col[u-1],ch_text[u-1],legend[u-1],u-1,aggregate_func[u-1],plots[u-1],
       plot_spacing[u-1],plot_height[u-1],v_align,tag_filter[u-1],
       plot_ticks[u-1],plot_color[u-1],plot_y_domain[u-1],
       ref_line[u-1],refl_pos[u-1],refl_color[u-1],refl_stroke_w[u-1],
       tagColor[u-1],renderHeat[u-1],text_font_size[u-1],
       label_font,label_angle,plot_filter[u-1],div_id,loci_links,ploidy_n,
       scatter_color_map_fn,show_links,seg_anno,directed_edges,
       links_colors,links_lg_x,links_lg_y,links_color_maps,win_scale,
       ann_h,display_chr[u-1],plot_shift[u-1],plot_legend_label[u-1],
       cat_legend_lab,plot_y_labs[u-1],plot_y_lab_x,plot_y_lab_y,
       plot_y_lab_size,interactivity,chr_text_x_shift);

u++;
}

/* end of rendering chromomap  */



} /* end of chromomap function */

function renderChromoMap(chData,chr_spacing,ploidy,ch_width,
        svg,nLoci,cnt,ch_gap,chr_color,heatmap,left_margin,lg_x,
        lg_y,heat_scale,labels,width,height,rng,heat_col,an_col,
        ch_text,legend,times,aggregate_func,
        plots,plot_spacing,plot_height,v,tag_f,plot_ticks,plot_color,
        plot_y_domain,ref_line,refl_pos,refl_color,refl_stroke_w,
        tagColor,renderHeat,text_font_size,label_font,label_angle,plot_f,
        div_id,loci_links,ploidy_n,scatter_color_map_fn,
        show_links,seg_anno,directed_edges,
        links_colors,links_lg_x,links_lg_y,links_color_maps,win_scale,
        ann_h,display_chr,plot_shift,plot_legend_label,
        cat_legend_lab,plot_y_labs,plot_y_lab_x,plot_y_lab_y,
        plot_y_lab_size,interactivity,chr_text_x_shift){


/* chromoMap render code */
if(cnt){

  // with centromere 

var allDatap =[];
var allDataq=[];
var allData = [];
for(i=0;i < nLoci.length;i++){

allDatap[i]=d3.range(nLoci[i].p);
allDataq[i]=d3.range(nLoci[i].q);
allData[i]=d3.range(nLoci[i].p+nLoci[i].q)
}

for(i=0;i< nLoci.length;i++){


if(display_chr){
/*  first cap */
posx=left_margin;
posy=y_val + plot_spacing + i*chr_spacing  + (ploidy - 1)*(ch_width +3 + plot_spacing) ;

if(!v){
svg.append("g")
.attr("transform", "translate("+posx+","+posy+")")
.append("path")
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",div_id+"-"+nLoci[i].name+"-"+1+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, true,false,true,false));

} else {

svg.append("g")
.attr("transform", "translate("+posx+","+posy+")")
.append("path")
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",div_id+"-"+nLoci[i].name+"-"+1+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, true,false,true,false));
}


/*  p arm */
if(!v){
svg.selectAll(".rects")
.data(allDatap[i].slice(1,allDatap[i].length-1))
.enter()
.append("rect")
.attr("y", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ))
.attr("height", ch_width)
.attr("x", function(d,i){return left_margin + loci_width + i*loci_width;})
.attr("width", loci_width)
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",function(d,j){ return div_id+"-"+nLoci[i].name+"-"+(j+2)+"-"+ploidy;});
} else {

svg.selectAll(".rects")
.data(allDatap[i].slice(1,allDatap[i].length-1))
.enter()
.append("rect")
.attr("x", y_val + i*chr_spacing + (ploidy - 1)*(ch_width +1 ))
.attr("width", ch_width)
.attr("y", function(d,i){return left_margin + i*loci_width;})
.attr("height", loci_width)
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",function(d,j){ return div_id+"-"+nLoci[i].name+"-"+(j+2)+"-"+ploidy;});


}


//labellings

if(!v){
svg.selectAll(".texts")
.data(allDatap[i])
.enter()
.append("text")
.attr("class",div_id+"-"+"labels")
.attr("y", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing )-2)
.attr("x", function(d,i){return left_margin + i*loci_width + (loci_width/2);})
.attr("id",function(d,j){ return div_id+"-"+"L"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;}).text("")
.attr("font-family", "sans-serif")
.attr("font-size", label_font)
.style("text-anchor", "start")
.attr("fill", "black");
} else {
svg.selectAll(".texts")
.data(allDatap[i])
.enter()
.append("text")
.attr("class","labels")
.attr("id",function(d,j) {return "L"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;}).text("")
.attr("font-family", "sans-serif")
.attr("font-size", "9px")
.style("text-anchor", "middle")
.attr("fill", "black").attr("transform",function(d,j){
return "translate(" + ((y_val + i*chr_spacing + (ploidy - 1)*(ch_width +1 ))-2) + ", "+(left_margin + j*loci_width)+") " + "rotate(-90)"});
}


/*  centromere */
posx=left_margin+ loci_width +(nLoci[i].p-2)*loci_width
//posy=y_val + plot_spacing + i*chr_spacing+  (ploidy - 1)*(ch_width +1 )
posy= y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing )

if(!v){
svg.append("g")
.attr("transform", "translate("+posx+","+posy+")")
.append("path")
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",div_id+"-"+nLoci[i].name+"-"+allDatap[i].length+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, false,true,false,true));
} else {
svg.append("g")
.attr("transform", "translate("+posy+","+posx+")")
.append("path")
.attr("fill", chr_color)
.attr("class","chLoc")
.attr("id",nLoci[i].name+"-"+allDatap[i].length+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, false,true,false,true));
}

posx=left_margin+ loci_width*2 +(nLoci[i].p-2)*loci_width
posy= y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing )

if(!v){
svg.append("g")
.attr("transform", "translate("+posx+","+posy+")")
.append("path")
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",div_id+"-"+nLoci[i].name+"-"+(nLoci[i].p+1)+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, true,false,true,false));
} else {
svg.append("g")
.attr("transform", "translate("+posy+","+posx+")")
.append("path")
.attr("fill", chr_color)
.attr("class","chLoc")
.attr("id",nLoci[i].name+"-"+(nLoci[i].p+1)+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, true,false,true,false));
}

/* short arm  */
if(!v){
svg.selectAll(".rects")
.data(allDataq[i].slice(1,allDataq[i].length-1))
.enter()
.append("rect")
.attr("y", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ))
.attr("height", ch_width)
.attr("x", function(d,j){return left_margin+ loci_width*3 + j*loci_width + (nLoci[i].p-2)*loci_width })
.attr("width", loci_width)
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",function(d,j){ return div_id+"-"+nLoci[i].name+"-"+(j+1+nLoci[i].p+1)+"-"+ploidy;});
} else {
svg.selectAll(".rects")
.data(allDataq[i].slice(1,allDataq[i].length-1))
.enter()
.append("rect")
.attr("x", y_val + i*chr_spacing+ (ploidy - 1)*(ch_width +1 ))
.attr("width", ch_width)
.attr("y", function(d,j){return left_margin+ loci_width*3 + j*loci_width + (nLoci[i].p-2)*loci_width })
.attr("height", loci_width)
.attr("fill", chr_color)
.attr("class","chLoc")
.attr("id",function(d,j){ return nLoci[i].name+"-"+(j+1+nLoci[i].p+1)+"-"+ploidy;});
}

//labellings

if(!v){
svg.selectAll(".texts")
.data(allDataq[i])
.enter()
.append("text")
.attr("class",div_id+"-"+"labels")
.attr("y", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing )-2)
.attr("x", function(d,j){return left_margin + j*loci_width + (nLoci[i].p)*loci_width + (loci_width/2);})
.attr("id",function(d,j){ return div_id+"-"+"L"+nLoci[i].name+"-"+(j+1+nLoci[i].p)+"-"+ploidy;}).text("")
.attr("font-family", "sans-serif")
.attr("font-size", label_font)
.style("text-anchor", "start")
.attr("fill", "black");
} else {
svg.selectAll(".texts")
.data(allDataq[i])
.enter()
.append("text")
.attr("class","labels")
.attr("id",function(d,j){return "L"+nLoci[i].name+"-"+(j+1+nLoci[i].p)+"-"+ploidy;}).text("")
.attr("font-family", "sans-serif")
.attr("font-size", "9px")
.style("text-anchor", "middle")
.attr("fill", "black").attr("transform",function(d,j){
return "translate(" + ((y_val + i*chr_spacing + (ploidy - 1)*(ch_width +1 ))-2) + ", "+(left_margin + j*loci_width+ (nLoci[i].p)*loci_width +2*arc_radius)+") " + "rotate(-90)"});
}





/*   final cap */


posx=left_margin+(nLoci[i].p-2)*loci_width + 3*loci_width + (nLoci[i].q-2)*loci_width;
posy=y_val + plot_spacing + i*chr_spacing+  (ploidy - 1)*(ch_width +3+ plot_spacing );
if(!v){
svg.append("g")
.attr("transform", "translate("+posx+","+posy+")")
.append("path")
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",div_id+"-"+nLoci[i].name+"-"+(nLoci[i].p +nLoci[i].q) +"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, false,true,false,true));
} else {
svg.append("g")
.attr("transform", "translate("+posy+","+posx+")")
.append("path")
.attr("fill", chr_color)
.attr("class","chLoc")
.attr("id",nLoci[i].name+"-"+(nLoci[i].p +nLoci[i].q) +"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, false,true,false,true));
}

if(!v){
if(ch_text){
 svg.append("text")
 .attr("x",left_margin*0.5)
 .attr("y",10 + y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ))
 .attr("font-family", "sans-serif")
 .attr("font-size", text_font_size)
 .attr("fill", "black")
 .style("text-anchor", "middle")
 .text(nLoci[i].name);
}
} else {
if(ch_text){
svg.append("text")
.attr("y",left_margin*0.75)
.attr("x",y_val + i*chr_spacing + (ploidy - 1)*(ch_width +2 ))
.attr("font-family", "sans-serif")
.attr("font-size", "9px")
.attr("fill", "black")
.style("text-anchor", "middle")
.text(nLoci[i].name);
   }
}


// end od display
}

if(plots=="bar"){
/* logic for bar*/
if(!v){

if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}

if(make_plot_y_labs){
svg.append("text")
 .attr("x",plot_y_lab_x)
 .attr("y",(y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g) + plot_y_lab_y )
 .attr("font-family", "sans-serif")
 .attr("font-size", plot_y_lab_size)
 .attr("fill", "black")
 .style("text-anchor", "middle")
 .text(plot_y_labs)
 .attr("class","plotylabs");

 svg.selectAll(".plotylabs")
 .attr("transform", function (d) {
  var xRot = d3.select(this).attr("x");
  var yRot = d3.select(this).attr("y");
  return `rotate(-90, ${xRot},  ${yRot} )`
  });
}

svg.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin - left_margin*0.5)+","+(y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g)+")";})
.attr("id",div_id+"-"+"bar-axis-"+nLoci[i].name+"-"+ploidy)
.attr("class",div_id+"-"+"plotxis");

svg.selectAll(".rects")
.data(allData[i])
.enter()
.append("rect")
.attr("y", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g)
.attr("height",plot_height)
.attr("x", function(d,i){return (left_margin + i*loci_width);})
.attr("width", loci_width)
.attr("fill", "white")
.attr("fill-opacity",0)
.attr("class",div_id+"-"+"barplot")
.attr("id",function(d,j){ return "bar-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;});
} else {
svg.selectAll(".rects")
.data(allData[i])
.enter()
.append("rect")
.attr("x", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing)
.attr("width",plot_height)
.attr("y", function(d,i){return (left_margin + i*loci_width);})
.attr("height", loci_width)
.attr("fill", "white")
.attr("class",div_id+"-"+"barplot")
.attr("id",function(d,j){ return div_id+"-"+"bar-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;});
}
}

/* logic for bar ends */
/* logic for scatter*/
if(plots=="scatter"){
if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}

sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +7+ plot_spacing ) - plot_spacing + plot_shift_g;

if(make_plot_y_labs){
svg.append("text")
.attr("x",plot_y_lab_x)
.attr("y",sc_y+plot_y_lab_y)
.attr("font-family", "sans-serif")
.attr("font-size", plot_y_lab_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_y_labs)
.attr("class","plotylabs");

svg.selectAll(".plotylabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});}



svg.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin - left_margin*0.5)+","+(sc_y)+")";})
.attr("id",div_id+"-"+"sc-axis-"+nLoci[i].name+"-"+ploidy)
.attr("class",div_id+"-"+"plotxis");
var test = svg.selectAll("scss")
.data(allData[i])
.enter()
.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin + j*loci_width)+","+sc_y+")";})
.attr("id",function(d,j){ return "sc-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;})
.attr("class",div_id+"-"+"scplot");


}

/* logic for scatter ends */
/* logic for 2d genome*/
if(plots=="2d"){
if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}

sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +7+ plot_spacing ) - plot_spacing + plot_shift_g;
if(make_plot_y_labs){
svg.append("text")
.attr("x",plot_y_lab_x)
.attr("y",sc_y+plot_y_lab_y)
.attr("font-family", "sans-serif")
.attr("font-size", plot_y_lab_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_y_labs)
.attr("class","plotylabs");

svg.selectAll(".plotylabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});}

svg.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin - left_margin*0.5)+","+(sc_y)+")";})
.attr("id",div_id+"-"+"g2d-axis-"+nLoci[i].name+"-"+ploidy)
.attr("class",div_id+"-"+"plotxis");
var test = svg.selectAll("g2dd")
.data(allData[i])
.enter()
.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin + j*loci_width)+","+sc_y+")";})
.attr("id",function(d,j){ return "g2d-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;})
.attr("class",div_id+"-"+"g2dplot");




}

/* logic for 2dg ends */
/* logic for tags*/
if(plots=="tags"){
if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}

sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g;

if(make_plot_y_labs){
svg.append("text")
.attr("x",plot_y_lab_x)
.attr("y",sc_y+plot_y_lab_y)
.attr("font-family", "sans-serif")
.attr("font-size", plot_y_lab_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_y_labs)
.attr("class","plotylabs");

svg.selectAll(".plotylabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});}
var test = svg.selectAll("tagss")
.data(allData[i])
.enter()
.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin + j*loci_width)+","+sc_y+")";})
.attr("id",function(d,j){ return "tags-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;})
.attr("class",div_id+"-"+"tags");


}
/* logic for tags ends */

/* logic for Lineplot*/
if(plots=="line"){
  if(!display_chr){
  plot_shift_g = ch_width + 10 + plot_shift;
  } else {
  plot_shift_g = 0;
  }
  
  sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g;
  
  if(make_plot_y_labs){
  svg.append("text")
  .attr("x",plot_y_lab_x)
  .attr("y",sc_y+plot_y_lab_y)
  .attr("font-family", "sans-serif")
  .attr("font-size", plot_y_lab_size)
  .attr("fill", "black")
  .style("text-anchor", "middle")
  .text(plot_y_labs)
  .attr("class","plotylabs");
  
  svg.selectAll(".plotylabs")
  .attr("transform", function (d) {
  var xRot = d3.select(this).attr("x");
  var yRot = d3.select(this).attr("y");
  return `rotate(-90, ${xRot},  ${yRot} )`
  });}
  var test = svg.selectAll("tagss")
  .data(allData[i])
  .enter()
  .append("g")
  .attr("transform", function(d,j){ return "translate("+(left_margin + j*loci_width)+","+sc_y+")";})
  .attr("id",function(d,j){ return "tags-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;})
  .attr("class",div_id+"-"+"tags");
  
  
  }
  /* logic for lineplot ends */

  

} //end of for

} else {
/**  without centromere*/
var allData =[];



for(i=0;i < nLoci.length;i++){

allData[i]=d3.range(nLoci[i].n);

}

//var deltaX2, deltaY2;
var dragHandler2 = d3.drag()
                    .on("start", 
                    function() {
                      d3.event.sourceEvent.stopPropagation();
                      //var current = d3.select(this);
                      //deltaX2 =  d3.event.dx;
                      //deltaY2 =  d3.event.dy;
                      d3.select(this).attr("stroke","red");
                    }
                    )
                    .on("drag", 
                    function(d) {
                      
                      d.x += d3.event.dx;
                      d.y += d3.event.dy;

                     d3.select(this).raise().attr("transform", function(d,i){
                           return "translate(" + [ d.x,d.y ] + ")"
                              });
                    }
                    )
                    .on("end", 
                    function() {
                      d3.select(this).attr("stroke",null);
                    }
                    );

  // adding the context menu
  
  

var clickFlag2=false;

for(i=0;i< nLoci.length;i++){

  
var labellings_per_ploidy = [false,true,true,true,true,true,true,true,true,true,true,true,true,true,true];
  

if(display_chr){

  if(labellings_per_ploidy[ploidy - 1]){
    plot_spacing = plot_spacing - 25;
  }

var chr_shift = (nLoci[i].n0*loci_width);
//left_margin = left_margin;
/*  first cap */
posx=left_margin + chr_shift;
posy=y_val + plot_spacing + i*chr_spacing  + (ploidy - 1)*(ch_width +3 + plot_spacing) ;


var chrGroup = svg.append("g")
                    .attr("class","chr_grp")
                    .datum({x:posx, y:posy})
                    .call(dragHandler2);

//var updateFill = function(){d3.select(this).style('fill',document.getElementById('fillcol001').value);}


chrGroup.on("contextmenu", function(d) {
  d3.event.stopPropagation();
  d3.event.preventDefault();
  var selctn = d3.select(this);
  var ctxtMenu =[];
ctxtMenu[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"ctxtMenu")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;	box-sizing: border-box; padding: 10px;	font: 12px sans-serif;		border: 0px;border-radius: 8px; ") ;
  const changeFill = function(){ selctn.selectAll("."+div_id+"-"+"chLoc").data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");}).attr("fill",document.getElementById('fillcol001').value);};
  const changeStroke = function(){ selctn.style('stroke',document.getElementById('strokecol001').value);};
  const changeFillChr = function(){ selctn.selectAll("."+div_id+"-"+"chLoc").attr("fill",document.getElementById('fillcol002').value);};
  //const closeThis = function(){ctxtMenu[ploidy-1].style("opacity", 0).style("z-index",-1);clickFlag2=false;}
  const closeThis = function(){ctxtMenu[ploidy-1].remove();clickFlag2=false;}
if(clickFlag2){
  ctxtMenu[ploidy-1].style("opacity", 0);
}else{
  ctxtMenu[ploidy-1].transition().duration(200)	.style("opacity", 1).style("z-index",1).style("background-color","#F5F5F5");	ctxtMenu[ploidy-1].style("height",100).style("width",300);
  ctxtMenu[ploidy-1].html("<div class='menuCtxt' style='display:flex;flex-direction:row;justify-content:center;align-items:center;padding:5px;gap:3px;'> <div class='menuItm' style='display:flex;flex-direction:column;justify-content:center;align-items:center;gap:3px;'><label> Fill(chr) </label> <input type='color' id='fillcol002' name='fillChr' value='#e66465' ><button class='fillChrButton' >update</button> </div><div class='menuItm' style='display:flex;flex-direction:column;justify-content:center;align-items:center;gap:3px;'><label> Fill(anno) </label> <input type='color' id='fillcol001' name='fillCol' value='#e66465'><button class='fillColButton' >update</button> </div><div class='menuItm' style='display:flex;flex-direction:column;justify-content:center;align-items:center;gap:3px;'> <label> Outline </label> <input type='color' id='strokecol001' name='strokeCol' value='#e66465'><button class='strokeColButton' >update</button></div><button class='closethis' style='border-radius:50%;'>x</button></div>").style("left", (d3.event.pageX) + "px")	.style("top", (d3.event.pageY - 18) + "px");
  document.getElementsByClassName("fillColButton")[0].addEventListener("click",changeFill);
  document.getElementsByClassName("strokeColButton")[0].addEventListener("click",changeStroke);
  document.getElementsByClassName("fillChrButton")[0].addEventListener("click",changeFillChr);
  document.getElementsByClassName("closethis")[0].addEventListener("click",closeThis);

} return clickFlag2 = !clickFlag2; });

if(!v){
chrGroup.append("g")
.append("path")
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",div_id+"-"+nLoci[i].name+"-"+1+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, true,false,true,false))
.transition()
.duration(500)
.attr("transform", "translate("+posx+","+posy+")");
} else {
svg.append("g")
.attr("transform", "translate("+posy+","+posx+")")
.append("path")
.attr("fill", chr_color)
.attr("class","chLoc")
.attr("id",nLoci[i].name+"-"+1+"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, true,false,true,false));
}

/*  p arm */

if(!v){
chrGroup.selectAll(".rects")
.data(allData[i].slice(1,allData[i].length-1))
.enter()
.append("rect")
.attr("height", ch_width)
.attr("width", loci_width)
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",function(d,j){ return div_id+"-"+nLoci[i].name+"-"+(j+2)+"-"+ploidy;})
.transition()
.duration(1000)
.attr("y", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ))
.attr("x", function(d,i){return (left_margin+ chr_shift + loci_width + i*loci_width);});

} else {
svg.selectAll(".rects")
.data(allData[i].slice(1,allData[i].length-1))
.enter()
.append("rect")
.attr("x", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ))
.attr("width", ch_width)
.attr("y", function(d,i){return (left_margin + loci_width + i*loci_width);})
.attr("height", loci_width)
.attr("fill", chr_color)
.attr("class","chLoc")
.attr("id",function(d,j){ return nLoci[i].name+"-"+(j+2)+"-"+ploidy;});
}



//labellings
if(!v){
svg.selectAll(".texts")
.data(allData[i])
.enter()
.append("text")
.attr("class",div_id+"-"+"labels")
.attr("y", (y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3 + plot_spacing))-2 )
.attr("x", function(d,i){return left_margin+ chr_shift + i*loci_width + (loci_width/2);})
.attr("id",function(d,j){return div_id+"-"+"L"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;}).text("")
.attr("font-family", "sans-serif")
.attr("font-size", label_font)
.style("text-anchor", "start")
.attr("fill", "black");
} else {
svg.selectAll(".texts")
.data(allData[i])
.enter()
.append("text")
.attr("class","labels")
.attr("id",function(d,j){return "L"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;}).text("")
.attr("font-family", "sans-serif")
.attr("font-size", "9px")
.style("text-anchor", "middle")
.attr("fill", "black").attr("transform",function(d,j){
return "translate(" + ((y_val + i*chr_spacing + (ploidy - 1)*(ch_width +1 ))-2) + ", "+(left_margin + j*loci_width)+") " + "rotate(-90)"});
}


/*   final cap */
posx=left_margin+ chr_shift+loci_width +(nLoci[i].n-2)*loci_width;
posy=y_val + plot_spacing + i*chr_spacing+  (ploidy - 1)*(ch_width +3+ plot_spacing );

if(!v){
chrGroup.append("g")
.append("path")
.attr("fill", chr_color)
.attr("class",div_id+"-"+"chLoc")
.attr("id",div_id+"-"+nLoci[i].name+"-"+ allData[i].length   +"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, false,true,false,true))
.transition()
.duration(500)
.attr("transform", "translate("+posx+","+posy+")");
} else {
svg.append("g")
.attr("transform", "translate("+posy+","+posx+")")
.append("path")
.attr("fill", chr_color)
.attr("class","chLoc")
.attr("id",nLoci[i].name+"-"+ allData[i].length   +"-"+ploidy)
.attr("d", rounded_rect(0, 0, loci_width, ch_width,ch_curve, false,true,false,true));
}

if(!v){

var chr_text_shift = [0,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80]
//var chr_text_x_shift = 65
if(ch_text){
svg.append("text")
.attr("x", left_margin*0.5  + chr_shift - chr_text_shift[ploidy - 1] + chr_text_x_shift) //left_margin*0.5
.attr("y",15+y_val + plot_spacing+ i*chr_spacing + (ploidy - 1)*(ch_width +3 + plot_spacing) )
.attr("font-family", "sans-serif")
.attr("font-size", text_font_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(nLoci[i].name);

}
} else {
if(ch_text){
 svg.append("text")
     .attr("y",left_margin*0.75)
     .attr("x",y_val + i*chr_spacing +arc_radius+ (ploidy - 1)*(ch_width +2 ))
     .attr("font-family", "sans-serif")
     .attr("font-size", "9px")
     .attr("fill", "black")
     .style("text-anchor", "middle")
     .text(nLoci[i].name);

}
}

/* end of visible chr*/
}

if(plots=="bar"){
/* logic for bar*/
if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}
if(!v){
if(make_plot_y_labs){
svg.append("text")
.attr("x",plot_y_lab_x)
.attr("y",(y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g) + plot_y_lab_y )
.attr("font-family", "sans-serif")
.attr("font-size", plot_y_lab_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_y_labs)
.attr("class","plotylabs");

svg.selectAll(".plotylabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});
}
svg.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin - left_margin*0.5)+","+(y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g )+")";})
.attr("id",div_id+"-"+"bar-axis-"+nLoci[i].name+"-"+ploidy)
.attr("class",div_id+"-"+"plotxis");

svg.selectAll(".rects")
.data(allData[i])
.enter()
.append("rect")
.attr("y", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g )
.attr("height",plot_height)
.attr("x", function(d,i){return (left_margin + i*loci_width);})
.attr("width", loci_width)
.attr("fill", "white")
.attr("fill-opacity",0)
.attr("class",div_id+"-"+"barplot")
.attr("id",function(d,j){ return "bar-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;});
} else {
svg.selectAll(".rects")
.data(allData[i])
.enter()
.append("rect")
.attr("x", y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing)
.attr("width",plot_height)
.attr("y", function(d,i){return (left_margin + i*loci_width);})
.attr("height", loci_width)
.attr("fill", "white")
.attr("class","barplot")
.attr("id",function(d,j){ return "bar-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;});
}
}

/* logic for bar ends */
/* logic for scatter*/
if(plots=="scatter"){
if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}

sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +7+ plot_spacing ) - plot_spacing + plot_shift_g;

if(make_plot_y_labs){
svg.append("text")
.attr("x",plot_y_lab_x)
.attr("y",sc_y+plot_y_lab_y)
.attr("font-family", "sans-serif")
.attr("font-size", plot_y_lab_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_y_labs)
.attr("class","plotylabs");

svg.selectAll(".plotylabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});}

svg.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin - left_margin*0.5)+","+(sc_y)+")";})
.attr("id",div_id+"-"+"sc-axis-"+nLoci[i].name+"-"+ploidy)
.attr("class",div_id+"-"+"plotxis");
var test = svg.selectAll("scss")
.data(allData[i])
.enter()
.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin + j*loci_width)+","+sc_y+")";})
.attr("id",function(d,j){ return "sc-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;})
.attr("class",div_id+"-"+"scplot");


}

/* logic for scatter ends */

/* logic for genome 2d*/
if(plots=="2d"){
if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}

sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +7+ plot_spacing ) - plot_spacing + plot_shift_g;

if(make_plot_y_labs){
svg.append("text")
.attr("x",plot_y_lab_x)
.attr("y",sc_y+plot_y_lab_y)
.attr("font-family", "sans-serif")
.attr("font-size", plot_y_lab_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_y_labs)
.attr("class","plotylabs");

svg.selectAll(".plotylabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});}

svg.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin - left_margin*0.5)+","+(sc_y)+")";})
.attr("id",div_id+"-"+"g2d-axis-"+nLoci[i].name+"-"+ploidy)
.attr("class",div_id+"-"+"plotxis");
var test = svg.selectAll("g2dd")
.data(allData[i])
.enter()
.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin + j*loci_width)+","+sc_y+")";})
.attr("id",function(d,j){ return "g2d-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;})
.attr("class",div_id+"-"+"g2dplot");


}

/* logic for genome 2d ends */
/* logic for tags*/
if(plots=="tags"){
if(!display_chr){
plot_shift_g = ch_width + 10 + plot_shift;
} else {
plot_shift_g = 0;
}

sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g;
if(make_plot_y_labs){
svg.append("text")
.attr("x",plot_y_lab_x)
.attr("y",sc_y+plot_y_lab_y)
.attr("font-family", "sans-serif")
.attr("font-size", plot_y_lab_size)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_y_labs)
.attr("class","plotylabs");

svg.selectAll(".plotylabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});}
var test = svg.selectAll("tagss")
.data(allData[i])
.enter()
.append("g")
.attr("transform", function(d,j){ return "translate("+(left_margin + j*loci_width)+","+sc_y+")";})
.attr("id",function(d,j){ return "tags-"+div_id+"-"+nLoci[i].name+"-"+(j+1)+"-"+ploidy;})
.attr("class",div_id+"-"+"tags");


}
/* logic for tage ends */

/* logic for lineplot*/
if(plots=="line"){
  if(!display_chr){
  plot_shift_g = ch_width + 10 + plot_shift;
  } else {
  plot_shift_g = 0;
  }
  
  sc_y=y_val + plot_spacing + i*chr_spacing + (ploidy - 1)*(ch_width +3+ plot_spacing ) - plot_spacing + plot_shift_g;
  if(make_plot_y_labs){
  svg.append("text")
  .attr("x",plot_y_lab_x)
  .attr("y",sc_y+plot_y_lab_y)
  .attr("font-family", "sans-serif")
  .attr("font-size", plot_y_lab_size)
  .attr("fill", "black")
  .style("text-anchor", "middle")
  .text(plot_y_labs)
  .attr("class","plotylabs");
  
  svg.selectAll(".plotylabs")
  .attr("transform", function (d) {
  var xRot = d3.select(this).attr("x");
  var yRot = d3.select(this).attr("y");
  return `rotate(-90, ${xRot},  ${yRot} )`
  });}


  var test = svg.append("g")
  .attr("transform", "translate("+(left_margin)+","+sc_y+")")
  .attr("id","line-"+div_id+"-"+nLoci[i].name+"-"+ploidy)
  .attr("class",div_id+"-"+"lines");
  
  
  }
  /* logic for lineplot ends */

} //end of for


} // end of render code

//var linking = false;
// code for linking loci
if(show_links){
if(ploidy==ploidy_n){
if(!seg_anno){
//console.log(loci_links[0].src_loci)
var loci_links_coords = [];
console.log(loci_links);
for(var tk=0;tk<loci_links.length;tk++){
let sx,sy,tx,ty;
if(d3.select("#"+loci_links[tk].src_loci).node().nodeName === "rect"){
sx = parseFloat(d3.select("#"+loci_links[tk].src_loci).attr("x"));
sy = parseFloat(d3.select("#"+loci_links[tk].src_loci).attr("y"));
} else if(d3.select("#"+loci_links[tk].src_loci).node().nodeName === "path"){
sx = (d3.select("#"+loci_links[tk].src_loci).node().getBoundingClientRect().x
- 10);
sy = (d3.select("#"+loci_links[tk].src_loci).node().getBoundingClientRect().y
- 10);
}
if(d3.select("#"+loci_links[tk].targ_loci).node().nodeName === "rect"){
tx = parseFloat(d3.select("#"+loci_links[tk].targ_loci).attr("x"));
ty = parseFloat(d3.select("#"+loci_links[tk].targ_loci).attr("y"));
} else if(d3.select("#"+loci_links[tk].targ_loci).node().nodeName === "path"){
tx = (d3.select("#"+loci_links[tk].targ_loci).node().getBoundingClientRect().x
- 10);
ty = (d3.select("#"+loci_links[tk].targ_loci).node().getBoundingClientRect().y
- 10);
}
//let sx_curve = d3.select("#chromap-chr2-1-1").node().getBoundingClientRect();
let gg = {};
if(sy > ty){
gg = {
x:[sx+loci_width/2,sy],
y:[tx+loci_width/2,ty+ch_width],
nm:loci_links[tk].lnk_nm,
dt:loci_links[tk].link_data
}
} else if(sy < ty){
gg = {
x:[sx+loci_width/2,sy+ch_width],
y:[tx+loci_width/2,ty],
nm:loci_links[tk].lnk_nm,
dt:loci_links[tk].link_data
}
} else {
gg = {
x:[sx+loci_width/2,sy],
y:[tx+loci_width/2,ty],
nm:loci_links[tk].lnk_nm,
dt:loci_links[tk].link_data
}
}

loci_links_coords.push(gg);
}
//console.log(loci_links_coords)

var dt_vals = [];
for(var i=0;i<loci_links_coords.length;i++){
dt_vals.push(loci_links_coords[i].dt)
}


//var directed_edges = false;
// directed edges
svg.append('defs').append('marker')
.attr('refX',13)
.attr('refY',0)
.attr('orient','auto')
.attr('markerWidth',10)
.attr('markerHeight',10)
.attr('xoverflow','visible')
.attr('viewBox','-0 -5 10 10')
.attr('id','arrowhead')
.append('svg:path')
.attr('d', 'M 0,-5 L 10 ,0 L 0,5')
.attr('fill', 'grey')
.style('stroke','none');

// string or number
const curve = d3.line().curve(d3.curveNatural);
var link_col_fn;
if(links_color_maps){
if(typeof dt_vals[0] == "string"){
var flags = [], uniq_vals = [];
for( var i=0; i<loci_links_coords.length; i++) {
if( flags[loci_links_coords[i].dt]) continue;
flags[loci_links_coords[i].dt] = true;
uniq_vals.push(loci_links_coords[i].dt);
}

uniq_vals = uniq_vals.sort();
link_col_fn = d3.scaleOrdinal()
.domain(uniq_vals)
.range(links_colors);

w22=uniq_vals.length*20;
// Create scale

var scalelinks = d3.scaleBand()
.domain(uniq_vals)
.range([2,w22]);

// Add scales to axis
var legendlinks_cat = d3.axisRight()
.scale(scalelinks).tickSizeOuter(0);


//Append group and insert axis
svg.append("g")
.attr("transform", "translate("+((width-35)-links_lg_x)+"," + ((height-100)-links_lg_y )+ ")")
.attr("class","scatter_legend")
.call(legendlinks_cat);

svg.select(".scatter_legend").style("font-size","15px");

rec_h22=w22/uniq_vals.length;
svg.selectAll(".rects")
.data(links_colors)
.enter()
.append("rect")
.attr("y", function(d,i){return (((height-99)-links_lg_y)+i*rec_h22);})
.attr("height", rec_h22)
.attr("x", function(d,i){return ((width-45)-links_lg_x);})
.attr("width", 10)
.attr("fill", function(d){return d;});




} else if(typeof dt_vals[0] == "number"){

link_col_fn = d3.scaleLinear()
.domain(d3.extent(dt_vals))
.range(links_colors);

link_legend_grad = svg.append("defs").append("linearGradient")
.attr("id", function(d){ return div_id+"-"+"linear-gradient-LINK-cont";})
.attr("x1", "0%")
.attr("y1", "0%")
.attr("x2", "0%")
.attr("y2", "100%")
.selectAll("stop")
.data( link_col_fn.range() )
.enter().append("stop")
.attr("offset", function(d,i) { return i/(link_col_fn.range().length-1); })
.attr("stop-color", function(d) { return d; });

svg.append('rect')
.attr("height",100).attr("width",10).attr("y",((height-99)-links_lg_y))
.attr("x",((width-45-(1*50))-links_lg_x))
.attr("class",function(d){return "lg_loci";})
.style("fill", "url(#"+div_id+"-"+"linear-gradient-LINK-cont)");

let locirng = [d3.min(dt_vals),d3.mean(dt_vals),d3.max(dt_vals)];
svg.selectAll(".texts")
.data(d3.range(3))
.enter()
.append("text")
.attr("class",div_id+"-"+"labels-links")
.text(function(i){ return locirng[i].toFixed(1);})
.attr("font-family", "sans-serif")
.attr("font-size", "9px")
.attr("fill", "black").attr("transform",function(d,j){
return "translate(" +((width-32-(1*50))-links_lg_x)+"," + ((height-94+(j*48))-links_lg_y  )+ ")"});

}
}
var links_tip  = d3.select("body").append("div")
.attr("class", div_id+"-"+"linkstooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 100px;height: 50px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;


svg.selectAll("paths")
.data(loci_links_coords)
.enter()
.append("path")
.attr("d", function(d){
if(d.x[1] != d.y[1]){
return d3.linkVertical()({
source: d.x,
target: d.y
});} else {
if(d.x[0]<d.y[0]){
let dist = d.y[0] - d.x[0];
let d_half = dist/2;
let d_quat = d_half/2;
let points = [];
points.push([d.x[0],d.x[1]]);
points.push([d.x[0]+d_quat,d.x[1]-ch_gap/2]);
points.push([d.x[0]+d_half,d.x[1]-ch_gap]);
points.push([d.x[0]+d_half+d_quat,d.x[1]-ch_gap/2]);
points.push([d.y[0],d.y[1]]);
return curve(points);
} else if(d.x[0]>d.y[0]){
let dist = d.x[0] - d.y[0];
let d_half = dist/2;
let d_quat = d_half/2;
let points = [];
points.push([d.x[0],d.x[1]]);
points.push([d.x[0]-d_quat,d.x[1]-ch_gap/2]);
points.push([d.x[0]-d_half,d.x[1]-ch_gap]);
points.push([d.x[0]-d_half-d_quat,d.x[1]-ch_gap/2]);
points.push([d.y[0],d.y[1]]);
return curve(points);

} else {
let points = [];
points.push([d.x[0],d.x[1]]);
points.push([d.x[0]-15,d.x[1]-ch_gap/2]);
points.push([d.x[0],d.x[1]-ch_gap]);
points.push([d.x[0]+15,d.x[1]-ch_gap/2]);
points.push([d.y[0],d.y[1]]);
return curve(points);
}
}
}).classed("link", true)
.attr("fill","none")
.attr("stroke-width","1.5")
.attr('marker-end',function(d){if(directed_edges){ return 'url(#arrowhead';} else { return "none"}})
.attr("stroke",function(d){   if(links_color_maps){ return link_col_fn(d.dt);}else{return links_colors}})
.on("mouseover", function(d) {
links_tip.transition().duration(200).style("opacity", .9).style("background-color","#F5F5F5").style("width",100+(d.dt.length*3)+'px');
links_tip.html("link: "+d.nm+"<br> value: "+d.dt).style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px"); d3.select(this).attr("stroke","black");
}).on("mouseout", function(d) {	links_tip.transition() .delay(1000).duration(500)	.style("opacity", 0);
d3.select(this).attr("stroke",function(d){if(links_color_maps){return link_col_fn(d.dt);}else{ return links_colors;}});});
}  else {

//logic for chords
//console.log(loci_links)
var loci_links_coords = [];

for(var tk=0;tk<loci_links.length;tk++){
let s1x,s1y,t1x,t1y,s2x,s2y,t2x,t2y;
if(d3.select("#"+loci_links[tk].src_loci).node().nodeName === "rect"){
s1x = parseFloat(d3.select("#"+loci_links[tk].src_loci).attr("x"));
s1y = parseFloat(d3.select("#"+loci_links[tk].src_loci).attr("y"));
} else if(d3.select("#"+loci_links[tk].src_loci).node().nodeName === "path"){
s1x = (d3.select("#"+loci_links[tk].src_loci).node().getBoundingClientRect().x
- 10);
s1y = (d3.select("#"+loci_links[tk].src_loci).node().getBoundingClientRect().y
- 10);
}
if(d3.select("#"+loci_links[tk].src_loci2).node().nodeName === "rect"){
s2x = parseFloat(d3.select("#"+loci_links[tk].src_loci2).attr("x"));
s2y = parseFloat(d3.select("#"+loci_links[tk].src_loci2).attr("y"));
} else if(d3.select("#"+loci_links[tk].src_loci2).node().nodeName === "path"){
s2x = (d3.select("#"+loci_links[tk].src_loci2).node().getBoundingClientRect().x
- 10);
s2y = (d3.select("#"+loci_links[tk].src_loci2).node().getBoundingClientRect().y
- 10);
}
if(d3.select("#"+loci_links[tk].targ_loci).node().nodeName === "rect"){
t1x = parseFloat(d3.select("#"+loci_links[tk].targ_loci).attr("x"));
t1y = parseFloat(d3.select("#"+loci_links[tk].targ_loci).attr("y"));
} else if(d3.select("#"+loci_links[tk].targ_loci).node().nodeName === "path"){
t1x = (d3.select("#"+loci_links[tk].targ_loci).node().getBoundingClientRect().x
- 10);
t1y = (d3.select("#"+loci_links[tk].targ_loci).node().getBoundingClientRect().y
- 10);
}
if(d3.select("#"+loci_links[tk].targ_loci2).node().nodeName === "rect"){
t2x = parseFloat(d3.select("#"+loci_links[tk].targ_loci2).attr("x"));
t2y = parseFloat(d3.select("#"+loci_links[tk].targ_loci2).attr("y"));
} else if(d3.select("#"+loci_links[tk].targ_loci2).node().nodeName === "path"){
t2x = (d3.select("#"+loci_links[tk].targ_loci2).node().getBoundingClientRect().x
- 10);
t2y = (d3.select("#"+loci_links[tk].targ_loci2).node().getBoundingClientRect().y
- 10);
}
//let sx_curve = d3.select("#chromap-chr2-1-1").node().getBoundingClientRect();
let gg = {};
if(s1y > t1y){
gg = {
s1:[s1x+loci_width/2,s1y],
t1:[t1x+loci_width/2,t1y+ch_width],
s2:[s2x+loci_width/2,s2y],
t2:[t2x+loci_width/2,t2y+ch_width],
nm:loci_links[tk].lnk_nm,
dt:loci_links[tk].link_data
}
} else if(s1y < t1y){
gg = {
s1:[s1x+loci_width/2,s1y+ch_width],
t1:[t1x+loci_width/2,t1y],
s2:[s2x+loci_width/2,s2y+ch_width],
t2:[t2x+loci_width/2,t2y],
nm:loci_links[tk].lnk_nm,
dt:loci_links[tk].link_data
}
} else {
gg = {
s1:[s1x+loci_width/2,s1y],
t1:[t1x+loci_width/2,t1y],
s2:[s2x+loci_width/2,s2y],
t2:[t2x+loci_width/2,t2y],
nm:loci_links[tk].lnk_nm,
dt:loci_links[tk].link_data
}
}

loci_links_coords.push(gg);
}
//console.log(loci_links_coords)

var dt_vals = [];
for(var i=0;i<loci_links_coords.length;i++){
dt_vals.push(loci_links_coords[i].dt)
}

// // string or number
const curve = d3.line().curve(d3.curveNatural);
var link_col_fn;
if(links_color_maps){
if(typeof dt_vals[0] == "string"){
var flags = [], uniq_vals = [];
for( var i=0; i<loci_links_coords.length; i++) {
if( flags[loci_links_coords[i].dt]) continue;
flags[loci_links_coords[i].dt] = true;
uniq_vals.push(loci_links_coords[i].dt);
}

link_col_fn = d3.scaleOrdinal()
.domain(uniq_vals)
.range(links_colors);

w22=uniq_vals.length*20;
// Create scale

var scalelinks = d3.scaleBand()
.domain(uniq_vals)
.range([2,w22]);

// Add scales to axis
var legendlinks_cat = d3.axisRight()
.scale(scalelinks).tickSizeOuter(0);


//Append group and insert axis
svg.append("g")
.attr("transform", "translate("+((width-35)-links_lg_x)+"," + ((height-100)-links_lg_y )+ ")")
.attr("class","scatter_legend")
.call(legendlinks_cat);

svg.select(".scatter_legend").style("font-size","15px");
rec_h22=w22/uniq_vals.length;
svg.selectAll(".rects")
.data(links_colors)
.enter()
.append("rect")
.attr("y", function(d,i){return (((height-99)-links_lg_y)+i*rec_h22);})
.attr("height", rec_h22)
.attr("x", function(d,i){return ((width-45)-links_lg_x);})
.attr("width", 10)
.attr("fill", function(d){return d;});


} else if(typeof dt_vals[0] == "number"){

link_col_fn = d3.scaleLinear()
.domain(d3.extent(dt_vals))
.range(links_colors);


link_legend_grad = svg.append("defs").append("linearGradient")
.attr("id", function(d){ return div_id+"-"+"linear-gradient-LINK-cont";})
.attr("x1", "0%")
.attr("y1", "0%")
.attr("x2", "0%")
.attr("y2", "100%")
.selectAll("stop")
.data( link_col_fn.range() )
.enter().append("stop")
.attr("offset", function(d,i) { return i/(link_col_fn.range().length-1); })
.attr("stop-color", function(d) { return d; });

svg.append('rect')
.attr("height",100).attr("width",10).attr("y",((height-99)-links_lg_y))
.attr("x",((width-45-(1*50))-links_lg_x))
.attr("class",function(d){return "lg_loci";})
.style("fill", "url(#"+div_id+"-"+"linear-gradient-LINK-cont)");

let locirng = [d3.min(dt_vals),d3.mean(dt_vals),d3.max(dt_vals)];
svg.selectAll(".texts")
.data(d3.range(3))
.enter()
.append("text")
.attr("class",div_id+"-"+"labels-links")
.text(function(i){ return locirng[i].toFixed(1);})
.attr("font-family", "sans-serif")
.attr("font-size", "9px")
.attr("fill", "black").attr("transform",function(d,j){
return "translate(" +((width-32-(1*50))-links_lg_x)+"," + ((height-94+(j*48))-links_lg_y  )+ ")"});

}
}
var links_tip  = d3.select("body").append("div")
.attr("class", div_id+"-"+"scattertooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 100px;height: 50px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;


svg.selectAll("paths")
.data(loci_links_coords)
.enter()
.append('path')
.attr("d", function(d){
if(d.s1[1] != d.t1[1]){
let chords = [];
chords.push(d3.linkVertical()({
source: d.s1,
target: d.t1
}));
chords.push(d3.linkVertical()({
source: d.t1,
target: d.t2
}));
chords.push(d3.linkVertical()({
source: d.t2,
target: d.s2
}));
//console.log(chords)
chords[0] = chords[0].split(",").splice(0,(chords[0].split(",").length-2)).join();
chords[0] = chords[0]+","
chords[1] = chords[1].split(",").splice(0,(chords[1].split(",").length - 2)).join();
chords[1] = chords[1].substring(1);
chords[1] = chords[1]+","
chords[2] = chords[2].substring(1);
//console.log(chords)
// chords.push(d3.linkVertical()({
//   source: d.s2,
//   target: d.s1
// }));
let full_chord_path = "";
for(var u=0;u<chords.length;u++){
full_chord_path += chords[u];
}
full_chord_path += "z";
console.log(full_chord_path);
return full_chord_path;

} else {
if(d.s1[0]<d.t1[0]){
let dist1 = d.t2[0] - d.s1[0];
let dist2 = d.t1[0] - d.s2[0];
let d1_half = dist1/2;
let d1_quat = d1_half/2;
let d2_half = dist2/2;
let d2_quat = d2_half/2;
let points1 = [];
let points2 = [];
points1.push([d.s1[0],d.s1[1]]);
points1.push([d.s1[0]+d1_quat,d.s1[1]-ch_gap/2]);
points1.push([d.s1[0]+d1_half,d.s1[1]-ch_gap]);
points1.push([d.s1[0]+d1_half+d1_quat,d.s1[1]-ch_gap/2]);
points1.push([d.t2[0],d.t2[1]]);
let out_cv = curve(points1);
points2.push([d.s2[0],d.s2[1]]);
points2.push([d.s2[0]+d2_quat,d.s2[1]-ch_gap/2]);
points2.push([d.s2[0]+d2_half,d.s2[1]-ch_gap + 5]);
points2.push([d.s2[0]+d2_half+d2_quat,d.s2[1]-ch_gap/2]);
points2.push([d.t1[0],d.t1[1]]);
let in_cv = curve(points2);
let lne1 = d3.line()([d.s1,d.s2]);
let lne2 = d3.line()([d.t1,d.t2]);
return lne1+in_cv+lne2+out_cv;
} else if(d.s1[0]>d.t1[0]){
let s_1 = d.s1;
let s_2 = d.s2;
let t_1 = d.t1;
let t_2 = d.t2;
let dist1 = s_2[0] - t_1[0];
let dist2 = s_1[0] - t_2[0];
let d1_half = dist1/2;
let d1_quat = d1_half/2;
let d2_half = dist2/2;
let d2_quat = d2_half/2;
let points1 = [];
let points2 = [];
points1.push([t_1[0],t_1[1]]);
points1.push([t_1[0]+d1_quat,t_1[1]-ch_gap/2]);
points1.push([t_1[0]+d1_half,t_1[1]-ch_gap]);
points1.push([t_1[0]+d1_half+d1_quat,t_1[1]-ch_gap/2]);
points1.push([s_2[0],s_2[1]]);
let out_cv = curve(points1);
points2.push([t_2[0],t_2[1]]);
points2.push([t_2[0]+d2_quat,t_2[1]-ch_gap/2]);
points2.push([t_2[0]+d2_half,t_2[1]-ch_gap + 5]);
points2.push([t_2[0]+d2_half+d2_quat,t_2[1]-ch_gap/2]);
points2.push([s_1[0],s_1[1]]);
let in_cv = curve(points2);
let lne1 = d3.line()([t_1,t_2]);
let lne2 = d3.line()([s_1,s_2])
return lne1+in_cv+lne2+out_cv;

} else {
let points = [];
points.push([d.s1[0],d.s1[1]]);
points.push([d.s1[0]-15,d.s1[1]-ch_gap/2]);
points.push([d.s1[0],d.s1[1]-ch_gap]);
points.push([d.s1[0]+15,d.s1[1]-ch_gap/2]);
points.push([d.t1[0],d.t1[1]]);
return curve(points);
}
}
}).classed("link", true)
.attr("fill",function(d){
if(d.s1[1] != d.t1[1]){
return link_col_fn(d.dt);
} else { return "none"}
})
.attr("stroke",function(d){  if(links_color_maps){ return link_col_fn(d.dt);}else { return links_colors;}})
.style("opacity",0.8)
.on("mouseover", function(d) {
links_tip.transition().duration(200).style("opacity", .9).style("background-color","#F5F5F5").style("width",100+(d.dt.length*3)+'px');
links_tip.html("link: "+d.nm+"<br> value: "+d.dt).style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px"); d3.select(this).attr("stroke","black");
}).on("mouseout", function(d) {	links_tip.transition() .delay(1000).duration(500)	.style("opacity", 0);
d3.select(this).attr("stroke",function(d){ if(links_color_maps){return link_col_fn(d.dt);}else{ return links_colors;}});});

//logic for chords ends

}
}}
// code ends fo linking loci

if(!heatmap){
//data reduction for anntation
var chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { var label =''; var tot = '<font color="grey"> Count: </font> '+v.length+"  </p></div><hr><div  style='float:top;height:80%;overflow-wrap: break-word;' > <font size='1'>";for(var i=0; i< v.length;i++){ label = label +"<a href="+v[i].hlink+" style='text-decoration:none;cursor:pointer;color:grey;pointer-events: all'>"+ v[i].name+"</a>" +" ,";    }; return tot+label.slice(0,-1);})
.entries(chData);

var chDataReducedCount = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  return v.length;})
.entries(chData);

print("HERRRRRRRR")
var chDataReducedOverlap = d3.nest()
.key(function(d) { return d.chrom; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].ch_start);};return bar;})
.entries(chData);

print(chDataReducedOverlap)
//print(chDataReducedCount);

var chDataRange = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a,b;for(var i=0; i< v.length;i++){ a = v[i].loci_start; b=v[i].loci_end;    };return ""+a+"-"+b+"bp ";})
.entries(chData);

/* adding the chromosome tootltips    */
for(var i = 0;i < chDataReduced.length;i++) {


chDataReduced[i].range = chDataRange[i].value;
chDataReduced[i].len = chDataReducedCount[i].value;
}



//data ready for rendering
var tip =[];
tip[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"chtooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;	box-sizing: border-box; padding: 3px;	font: 12px sans-serif;		border: 0px;border-radius: 8px;pointer-events: none;") ;

var tip2 = [];
tip2[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"chtooltip2")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 180px;height: 140px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;
var clickFlag=false;
if(interactivity){
  
d3.selectAll("."+div_id+"-"+"chLoc").data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");}).on("mouseover", function(d) {tip[ploidy-1].transition().duration(200)	.style("opacity", .9).style("background-color","#F5F5F5");	tip[ploidy-1].style("height",100+d.len*1.5).style("width",190+d.len);
tip[ploidy-1].html("<div  style='float:top;height:10%;color:black;font-size:12px;'><p>"+d.range+"</p></div><hr><div  style='float:top;height:10%;'><p ><b>"+d.value+"</font></div>").style("left", (d3.event.pageX) + "px")	.style("top", (d3.event.pageY - 18) + "px");}).on("mouseout", function(d) {	tip[ploidy-1].transition() .delay(1000).duration(500)	.style("opacity", 0);}).attr("fill", an_col)
.on("click", function(d) {


if(clickFlag){
tip2[ploidy-1].style("opacity", 0);
}else{
tip2[ploidy-1].transition().duration(200)	.style("opacity", 1).style("background-color","#F5F5F5");	tip2[ploidy-1].style("height",100+d.len*1.5).style("width",190+d.len);
tip2[ploidy-1].html("<div  style='float:top;height:10%;color:black;font-size:12px;'><p>"+d.range+"</p></div><hr><div  style='float:top;height:10%;'><p ><b>"+d.value+"</font></div>").style("left", (d3.event.pageX) + "px")	.style("top", (d3.event.pageY - 18) + "px");
} return clickFlag = !clickFlag; });
} else {
// without interactivity
d3.selectAll("."+div_id+"-"+"chLoc")
.data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.attr("fill", an_col);

//end of without interactivity
}



if(labels){

var chDataLabels = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a=[];for(var i=0; i< v.length;i++){ a.push(v[i].name);   };return a[0];})
.entries(chData);
var chDataLabel = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a=[];for(var i=0; i< v.length;i++){ a.push(v[i].label);   };return a[0];})
.entries(chData);
for(var i = 0;i < chDataReduced.length;i++) {



chDataReduced[i].labels =chDataLabels[i].value;
chDataReduced[i].label_id =chDataLabel[i].value;

}

d3.selectAll("."+div_id+"-"+"labels").data(chDataReduced, function(d) {return (d && d.label_id)||d3.select(this).attr("id");})
.text(function(d){ return d.labels;})
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(${label_angle}, ${xRot},  ${yRot} )`
});
}




} else {

//heatmap code
if(heat_scale=="numeric"){


// if(times==1 || times==2){
//   rng=[2.5,3.8];
// }

colors[times] = d3.scaleLinear()
.domain(rng).nice()
.range(heat_col);

print("insideeee heat ")
print(rng)
print(heat_col)







// Create scale

/*    var scale2 = d3.scaleLinear()
.domain(rng).nice()
.range([2, 100]);

// Add scales to axis
legend2 = d3.axisRight()
.scale(scale2)
.ticks(3);
if(legend){
//Append group and insert axis
svg.append("g")
.attr("transform", "translate("+((width-35-(times*50))-lg_x)+"," + ((height-100)-lg_y )+ ")")
.call(legend2);

}  */
//console.log(colors[times].range());
chLinGradV[times] = svg.append("defs").append("linearGradient")
        .attr("id", function(d){ return div_id+"-"+"linear-gradient-chromH-"+times;})
        .attr("x1", "0%")
        .attr("y1", "0%")
        .attr("x2", "0%")
        .attr("y2", "100%")
        .selectAll("stop")
.data( colors[times].range().reverse() )
.enter().append("stop")
.attr("offset", function(d,i) { return i/(colors[times].range().reverse().length-1); })
.attr("stop-color", function(d) { return d; });


if(legend){
svg.append('rect')
.attr("height",100).attr("width",12).attr("y",((height-99)-lg_y)).attr("x",((width-45+(times*100))-lg_x)).attr("class",function(d){return "lg"+times;})
.style("fill", "url(#"+div_id+"-"+"linear-gradient-chromH-"+times+")");

var rng2;
// if(rng.length==3){
//   rng2=rng;
// }else if (rng.length==2) {
//   rng2=[rng[0],d3.mean(rng),rng[1]]
// }

if(rng.length==3){
rng2=[rng[2],rng[1],rng[0]];
}else if (rng.length==2) {
rng2=[rng[1],d3.mean(rng),rng[0]]
}

svg.selectAll(".texts")
.data(d3.range(3))
.enter()
.append("text")
.attr("class",div_id+"-"+"labels")
.text(function(i){ return rng2[i].toFixed(1);})
.attr("font-family", "sans-serif")
.attr("font-size", "12px")
.attr("fill", "black").attr("transform",function(d,j){
return "translate(" +((width-32+(times*100))-lg_x + 2)+"," + ((height-90+(j*44))-lg_y  )+ ")"});


svg.append("text")
.attr("x",((width-45+(times*100))-lg_x - 10))
.attr("y",((height-99)-lg_y+50) )
.attr("font-family", "sans-serif")
.attr("font-size", 15)
.attr("fill", "black")
.style("text-anchor", "middle")
.text(plot_legend_label)
.attr("class","plotleglabs");

svg.selectAll(".plotleglabs")
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(-90, ${xRot},  ${yRot} )`
});

}

} else{ if(heat_scale=="categorical"){




colors[times] = d3.scaleOrdinal()
.domain(rng)
.range(heat_col);


w=rng.length*20;
// Create scale

var scale2 = d3.scaleBand()
.domain(rng)
.range([2,w]);

// Add scales to axis
var legend2 = d3.axisRight()
.scale(scale2).tickSizeOuter(0);

if(legend){
//Append group and insert axis
svg.append("g")
.attr("transform", "translate("+((width-35)-lg_x)+"," + ((height-100)-lg_y )+ ")")
.call(legend2).attr("id","glgnd");

svg.select("#glgnd").style("font-size","20px");

// LEGEND LABEL
//ading the label for the LEGEND
svg.append("text")
.attr("x",(width-35)-lg_x)
.attr("y",((height-100)-lg_y - 15 ) )
.attr("font-family", "sans-serif")
.attr("font-size", 20)
.attr("fill", "black")
.style("text-anchor", "start")
.text(cat_legend_lab)
.attr("class","sctrlbl2");



rec_h=w/rng.length;

svg.selectAll(".rects")
.data(heat_col)
.enter()
.append("rect")
.attr("y", function(d,i){return (((height-99)-lg_y)+i*rec_h);})
.attr("height", rec_h)
.attr("x", function(d,i){return ((width-45)-lg_x);})
.attr("width", 10)
.attr("fill", function(d){return d;});
}




}
}




/* creating final data    */



if(heat_scale=="numeric"){

/* aggregate functions*/
var tag;
if(aggregate_func=="avg"){
chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.mean(v, function(d) { return d.data; });})
.entries(chData);

tag="avg";
} else if(aggregate_func=="sum"){

chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.sum(v, function(d) { return d.data; });})
.entries(chData);
tag="sum";

} else if(aggregate_func=="min"){

chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.min(v, function(d) { return d.data; });})
.entries(chData);
tag="min";

} else if(aggregate_func=="max"){

chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.max(v, function(d) { return d.data; });})
.entries(chData);
tag="max";

} else if(aggregate_func=="count"){

chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return v.length;})
.entries(chData);
tag="count";

} else if (aggregate_func=="median"){
chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.median(v, function(d) { return d.data; });})
.entries(chData);
tag="median";
} else if (aggregate_func=="none"){
chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.mean(v, function(d) { return d.data; });})
.entries(chData);
}
/* aggregate functions*/





var chDataReducedMin = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.min(v, function(d) { return d.data; });})
.entries(chData);
var chDataReducedMax = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return d3.max(v, function(d) { return d.data; });})
.entries(chData);
var chDataReducedBarData = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].data);};return bar;})
.entries(chData);
var chDataReducedBarLabel = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].name);};return bar;})
.entries(chData);

var chDataReducedCategory = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].cate);};return bar;})
.entries(chData);

var chDataReducedBarLink = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].hlink);};return bar;})
.entries(chData);
var chDataReducedCount = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return v.length;})
.entries(chData);
var chDataRange = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a,b;for(var i=0; i< v.length;i++){ a = v[i].loci_start; b=v[i].loci_end;    };return ""+a+"-"+b+"bp ";})
.entries(chData);

var chDataReducedForStatic = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { var label =''; var tot = '<font color="grey"> Count: </font> '+v.length+"  </p></div><hr><div  style='float:top;height:80%;overflow-wrap: break-word;' > <font size='1'>";for(var i=0; i< v.length;i++){ label = label +"<a href="+v[i].hlink+" style='text-decoration:none;cursor:pointer;color:grey;pointer-events: all'>"+ v[i].name+"</a>" +" ,";    }; return tot+label.slice(0,-1);})
.entries(chData);
for(var i = 0;i < chDataReduced.length;i++) {


chDataReduced[i].range = chDataRange[i].value;
chDataReduced[i].min = chDataReducedMin[i].value;
chDataReduced[i].max = chDataReducedMax[i].value;
chDataReduced[i].bar = chDataReducedBarData[i].value;
chDataReduced[i].label = chDataReducedBarLabel[i].value;
chDataReduced[i].count = chDataReducedCount[i].value;
chDataReduced[i].hlink = chDataReducedBarLink[i].value;
chDataReduced[i].static = chDataReducedForStatic[i].value;
chDataReduced[i].category = chDataReducedCategory[i].value;
}

if(plots == "line"){

  var chDataReducedLineLoci = d3.nest()
.key(function(d) { return d.chr_id; })
.rollup(function(v){  arr = []; for(var i=0; i< v.length;i++) {  arr.push({loci:v[i].loci,data:v[i].data})  }; 
                 return d3.nest().key(function(d) {return d.loci;}).rollup(function(v) { return d3.mean(v, function(d) { return d.data; });}).entries(arr);  })
.entries(chData);

for (var i = 0; i < chDataReducedLineLoci.length; i++) {

  let l = chDataReducedLineLoci[i].value;

  for(var j = 0; j < l.length;j++){
    l[j].win = parseFloat(l[j].key.split("-")[2]);

  }
  }

}

print(chDataReducedLineLoci)


delete chDataRange;
delete chDataReducedMin;
delete chDataReducedMax;
delete chDataReducedBarData;
delete chDataReducedBarLabel;
delete chDataReducedBarLink;
delete chDataReducedForStatic;
delete chDataReducedCategory;

if(plots=="bar"){
/*barplot*/
var BarData=JSON.parse(JSON.stringify(chDataReduced));
for (var i = 0; i < BarData.length; i++) {

BarData[i].key = "bar-"+BarData[i].key;
}


// Add Y axis
if(parseFloat(plot_y_domain[0]) === 0 && parseFloat(plot_y_domain[1]) === 0){
var y = d3.scaleLinear()
.domain(d3.extent(BarData, function(d) { return d.value} ))
.range([ plot_height,0]).nice();
} else {
var y = d3.scaleLinear()
.domain(plot_y_domain)
.range([ plot_height,0]).nice();
}
for(i=0;i< nLoci.length;i++){
//console.log("axis-"+nLoci[i].name+"-"+ploidy)
svg.select("#"+div_id+"-"+"bar-axis-"+nLoci[i].name+"-"+ploidy).call(d3.axisLeft().scale(y).ticks(plot_ticks));


}

//console.log(y(1));
}
/*barplot*/
//console.log(BarData)
/*scatterplot*/
if(plots=="scatter"){
var ScData=JSON.parse(JSON.stringify(chDataReduced));
for (var i = 0; i < ScData.length; i++) {

ScData[i].key = "sc-"+ScData[i].key;
}

//console.log(plot_y_domain);
// Add Y axis
if(parseFloat(plot_y_domain[0]) === 0 && parseFloat(plot_y_domain[1]) === 0){
var scy = d3.scaleLinear()
.domain([d3.min(ScData, function(d) { return d.min} ),d3.max(ScData, function(d) { return d.max} )])
.range([ (plot_height - 5),0]).nice();
} else {
var scy = d3.scaleLinear()
.domain(plot_y_domain)
.range([ (plot_height - 5),0]).nice();
}
for(i=0;i< nLoci.length;i++){

svg.select("#"+div_id+"-"+"sc-axis-"+nLoci[i].name+"-"+ploidy).call(d3.axisLeft().scale(scy).ticks(plot_ticks))
}

}

/*2d genome*/
if(plots=="2d"){
var g2dData=JSON.parse(JSON.stringify(chDataReduced));
for (var i = 0; i < g2dData.length; i++) {

g2dData[i].key = "g2d-"+g2dData[i].key;
}

//console.log(plot_y_domain);
// Add Y axis
if(parseFloat(plot_y_domain[0]) === 0 && parseFloat(plot_y_domain[1]) === 0){
var g2dy = d3.scaleLinear()
.domain([d3.min(g2dData, function(d) { return d.min} ),d3.max(g2dData, function(d) { return d.max} )])
.range([ (plot_height - 5 ),0]).nice();
} else {
var g2dy = d3.scaleLinear()
.domain(plot_y_domain)
.range([ (plot_height - 5 ),0]).nice();
}
for(i=0;i< nLoci.length;i++){

svg.select("#"+div_id+"-"+"g2d-axis-"+nLoci[i].name+"-"+ploidy).call(d3.axisLeft().scale(g2dy).ticks(plot_ticks))
}

}

// adding logic for the axis.


/*scatterplot*/
if(plots=="tags"){
var tagDatatemp=JSON.parse(JSON.stringify(chDataReduced));
for (var i = 0; i < tagDatatemp.length; i++) {

tagDatatemp[i].key = "tags-"+tagDatatemp[i].key;
}

switch(tag_f[0]){
case "eq":
var tagData = tagDatatemp.filter(function(d){ return d.value == tag_f[1] });
break;
case "gt":
var tagData = tagDatatemp.filter(function(d){ return d.max > tag_f[1] });
break;
case "gte":
var tagData = tagDatatemp.filter(function(d){ return d.max >= tag_f[1] });
break;
case "lt":
var tagData = tagDatatemp.filter(function(d){ return d.min < tag_f[1] });
break;
case "lte":
var tagData = tagDatatemp.filter(function(d){ return d.min <= tag_f[1] });
break;
case "gtalt":
var tagData = tagDatatemp.filter(function(d){ return d.max > tag_f[1] && d.min < tag_f[2] });
break;
case "gtealte":
var tagData = tagDatatemp.filter(function(d){ return d.max >= tag_f[1] && d.min <= tag_f[2] });
break;
case "gtolt":
var tagData = tagDatatemp.filter(function(d){ return d.max > tag_f[1] || d.min < tag_f[2] });
break;
case "gteolte":
var tagData = tagDatatemp.filter(function(d){ return d.max >= tag_f[1] || d.min <= tag_f[2] });
break;
case "none":
var tagData=JSON.parse(JSON.stringify(tagDatatemp));
default:
console.log("nothing.")

}

//var tagData = tagDatatemp.filter(function(d){ return d.min > 55 })
//console.log(tagData.filter(function(d){ return d.min > 55 }))
}

/* */
var tip = [];
tip[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"chtooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 180px;height: 140px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;

var tip2 = [];
tip2[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"chtooltip2")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 180px;height: 140px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;
var clickFlag=false;


if(renderHeat){

if(interactivity) {
d3.selectAll("."+div_id+"-"+"chLoc").data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");}).on("mouseover", function(d) {tip[ploidy-1].transition().duration(200)	.style("opacity", .9).style("background-color","#F5F5F5").style("width",180+d.count*(d3.max(d.label).length*3)+'px');
tip[ploidy-1].html("<div  style='float:top;position:relative;height:35%;' ><div style='border-radius: 5px;float:left;position:relative;height:55px;width:30%;background-color:"+colors[times](d.value)+
";' ></div><div  style='float:left;position:relative;width:70%;' ><p ><font size='1' color='grey' > "+tag+": </font> <font size='1'  ><b> "+(d.value).toFixed(2)+
" <br><font size='1' color='grey' > Min.: </font>"+(d.min).toFixed(2)+" <br><font size='1' color='grey' > Max.: </font> "+(d.max).toFixed(2)+
"</font><br><font size='1' color='grey' > Count: </font><font size='1'>"+d.count+"</p></font></div></div><div id ='"+div_id+"-microBar"+ploidy+"' style='float:top;position:relative;height:35%'></div><br><div  style='float:top;height:20%;position:relative;color:'black"+
";'><font size='2' >"+d.range+"</font></p></div>")
.style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px"); var svgWidth = 175+d.count*(d3.max(d.label).length)*3;var svgHeight = 50;var barPadding = 5;var barWidth = (svgWidth / d.bar.length);var mysvg = d3.select('#'+div_id+'-microBar'+ploidy).append("svg").attr("width",svgWidth).attr("heigth",svgHeight);var barChart = mysvg.selectAll('rect').data(d.bar).enter().append('rect').attr('height',svgHeight-40).attr('width', barWidth).attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';}).style("fill",function(d){return colors[times](d);});
t_data=[];
for(j=0;j<d.label.length;j++){
t_data.push({lb:d.label[j],hl:d.hlink[j]}); }
var ttext=mysvg.selectAll(".node").data(t_data).enter().append("svg:a").attr("xlink:href", function(d,i){ return d.hl;})
.append("svg:text").attr("y",svgHeight*0.4).attr("transform","rotate(90)")
.attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';})
.style("pointer-events", "all").style("cursor","pointer").attr("font-size", "9px")
.text(function(d,i){return d.lb;}); })
.on("mouseout", function(d) {	tip[ploidy-1].transition() .delay(1000).duration(500)	.style("opacity", 0);}).style("visibility","visible").style("fill",function(d){ return colors[times](d.value);})
.on("click", function(d) {


if(clickFlag){
tip2[ploidy-1].style("opacity", 0);
}else{
tip2[ploidy-1].style("z-index",99999999)	.style("opacity", 1).style("background-color","#F5F5F5").style("width",180+d.count*(d3.max(d.label).length*3)+'px');
tip2[ploidy-1].html("<div  style='float:top;position:relative;height:35%;' ><div style='border-radius: 5px;float:left;position:relative;height:55px;width:30%;background-color:"+colors[times](d.value)+
";' ></div><div  style='float:left;position:relative;width:70%;' ><p ><font size='1' color='grey' > "+tag+": </font> <font size='1'  ><b> "+(d.value).toFixed(2)+
" <br><font size='1' color='grey' > Min.: </font>"+(d.min).toFixed(2)+" <br><font size='1' color='grey' > Max.: </font> "+(d.max).toFixed(2)+
"</font><br><font size='1' color='grey' > Count: </font><font size='1'>"+d.count+"</p></font></div></div><div id ='"+div_id+"-microBar2"+ploidy+"' style='float:top;position:relative;height:35%'></div><br><div  style='float:top;height:20%;position:relative;color:'black"+
";'><font size='2' >"+d.range+"</font></p></div>")
.style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px"); var svgWidth = 175+d.count*(d3.max(d.label).length)*3;var svgHeight = 50;var barPadding = 5;var barWidth = (svgWidth / d.bar.length);var mysvg2 = d3.select('#'+div_id+'-microBar2'+ploidy).append("svg").attr("width",svgWidth).attr("heigth",svgHeight);var barChart = mysvg2.selectAll('rect').data(d.bar).enter().append('rect').attr('height',svgHeight-40).attr('width', barWidth).attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';}).style("fill",function(d){return colors[times](d);});
t_data=[];
for(j=0;j<d.label.length;j++){
t_data.push({lb:d.label[j],hl:d.hlink[j]}); }
var ttext=mysvg2.selectAll(".node").data(t_data).enter().append("svg:a").attr("xlink:href", function(d,i){ return d.hl;}).append("svg:text").attr("y",svgHeight*0.4).attr("transform","rotate(90)").attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';}).style("pointer-events", "all").style("cursor","pointer").attr("font-size", "9px").text(function(d,i){return d.lb;});
} return clickFlag = !clickFlag; });
} else {
//without interacitivuty
d3.selectAll("."+div_id+"-"+"chLoc")
.data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.style("visibility","visible")
.style("fill",function(d){ return colors[times](d.value);});


//without interactivity end
}


/*end of renderHeat*/ } else {

if(interactivity){
d3.selectAll("."+div_id+"-"+"chLoc").data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");}).on("mouseover", function(d) {tip[ploidy-1].transition().duration(200)	.style("opacity", .9).style("background-color","#F5F5F5");	tip[ploidy-1].style("height",100+d.len*1.5).style("width",190+d.len);
tip[ploidy-1].html("<div  style='float:top;height:10%;color:black;font-size:12px;'><p>"+d.range+"</p></div><hr><div  style='float:top;height:10%;'><p ><b>"+d.static+"</font></div>").style("left", (d3.event.pageX) + "px")	.style("top", (d3.event.pageY - 18) + "px");}).on("mouseout", function(d) {	tip[ploidy-1].transition() .delay(1000).duration(500)	.style("opacity", 0);}).attr("fill", an_col)
.on("click", function(d) {


if(clickFlag){
tip2[ploidy-1].style("opacity", 0);
}else{
tip2[ploidy-1].transition().duration(200)	.style("opacity", 1).style("background-color","#F5F5F5");	tip2[ploidy-1].style("height",100+d.len*1.5).style("width",190+d.len);
tip2[ploidy-1].html("<div  style='float:top;height:10%;color:black;font-size:12px;'><p>"+d.range+"</p></div><hr><div  style='float:top;height:10%;'><p ><b>"+d.static+"</font></div>").style("left", (d3.event.pageX) + "px")	.style("top", (d3.event.pageY - 18) + "px");
} return clickFlag = !clickFlag; });
} else {
//without interactvity
d3.selectAll("."+div_id+"-"+"chLoc")
.data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.attr("fill", an_col);

//without interactivity ends
}

/*end of renderHeat else*/}

if(labels){

var chDataLabels = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a=[];for(var i=0; i< v.length;i++){ a.push(v[i].name);   };return a[0];})
.entries(chData);

var chDataLabel = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a=[];for(var i=0; i< v.length;i++){ a.push(v[i].label);   };return a[0];})
.entries(chData);
for(var i = 0;i < chDataReduced.length;i++) {


chDataReduced[i].label_id=chDataLabel[i].value;
chDataReduced[i].labels =chDataLabels[i].value;

}
//console.log(chDataReduced);
delete chDataLabels;

d3.selectAll("."+div_id+"-"+"labels").data(chDataReduced, function(d) {return (d && d.label_id)||d3.select(this).attr("id");})
.text(function(d){ return d.labels;})
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(${label_angle}, ${xRot},  ${yRot} )`
});
}


if(plots=="bar"){
/*bar*/

d3.selectAll("."+div_id+"-"+"barplot").data(BarData, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.attr("fill",function(d){
//if(d.value >= 0){ return plot_color;} else {return "red";}
switch(plot_f[0]){
case "eq":
if(d.value == plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gt":
if(d.value > plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gte":
if(d.value >= plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "lt":
if(d.value < plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "lte":
if(d.value <= plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gtalt":
if(d.value > plot_f[1] && d.value < plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gtealte":
if(d.value >= plot_f[1] && d.value <= plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gtolt":
if(d.value > plot_f[1] || d.value < plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gteolte":
if(d.value >= plot_f[1] || d.value <= plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "none":
return plot_color;
default:
return plot_color;

}

})
.attr("fill-opacity",1)
.attr("y", function(d){ return d3.select(this).attr('y') - 0.001+y(d.value);})
.attr("height",function(d) { return plot_height-y(d.value); });
if(ref_line){
for(i=0;i< nLoci.length;i++){

svg.select("#"+div_id+"-"+"bar-axis-"+nLoci[i].name+"-"+ploidy).append("line")
.attr("x1", 0)
.attr("x2", left_margin +  (win_scale*loci_width))
.attr("y1",  y(refl_pos))
.attr("y2", y(refl_pos))
.attr("stroke", refl_color)
.attr("stroke-dasharray","5,5")
.attr("stroke-width",refl_stroke_w)
.attr("class",div_id+"-"+"overlines");


}}
}
/**/
/*scatter*/
if(plots=="scatter"){
//console.log(ScData);
var sc_tip = [];
sc_tip[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"scattertooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 100px;height: 50px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;

d3.selectAll("."+div_id+"-"+"scplot").data(ScData, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.selectAll("dots")
.data(function(d,i){
var sc_data = [];
for(var l=0;l<d.bar.length;l++){
var objj = {"br":d.bar[l],"nm":d.label[l],"cte":d.category[l]};
sc_data.push(objj);
}
return sc_data;
})
.enter()
.append("circle")
.attr("cx", loci_width/2 )
.attr("cy", function (d) { return scy(d.br); } )
.attr("r", 1.5)
.style("fill", function(d){
//if(d.br >= 0){ return plot_color;} else {return "red";}
//var scatter_fill;
switch(plot_f[0]){
case "eq":
if(d.br == plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gt":
if(d.br > plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gte":
if(d.br >= plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "lt":
if(d.br < plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "lte":
if(d.br <= plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gtalt":
if(d.br > plot_f[1] && d.br < plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gtealte":
if(d.br >= plot_f[1] && d.br <= plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gtolt":
if(d.br > plot_f[1] || d.br < plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gteolte":
if(d.br >= plot_f[1] || d.br <= plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "none":
return plot_color;
case "col":
return scatter_color_map_fn(d.cte);
default:
return plot_color;

}

})
.on("mouseover", function(d) {
sc_tip[ploidy-1].transition().duration(200).style("opacity", .9).style("background-color","#F5F5F5").style("width",100+(d.nm.length*3)+'px');
sc_tip[ploidy-1].html("name: "+d.nm+"<br> value: "+d.br.toFixed(2)).style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px");
}).on("mouseout", function(d) {	sc_tip[ploidy-1].transition() .delay(1000).duration(500)	.style("opacity", 0);});
if(ref_line){
for(i=0;i< nLoci.length;i++){

svg.select("#"+div_id+"-"+"sc-axis-"+nLoci[i].name+"-"+ploidy).append("line")
.attr("x1", 0)
.attr("x2", left_margin +  (win_scale*loci_width))
.attr("y1",  scy(refl_pos))
.attr("y2", scy(refl_pos) )
.attr("stroke", refl_color)
.attr("stroke-dasharray","5,5")
.attr("stroke-width",refl_stroke_w)
.attr("class",div_id+"-"+"overlines");


}}
}
/**/

var gnbox = { 
draw: function(context, size){
let s = Math.sqrt(size)/2;
let h = ann_h;
context.moveTo(s,h);
context.lineTo(s,-h);
context.lineTo(-s,-h);
context.lineTo(-s,h);
context.closePath();
}
}

/*2d genome plot*/
if(plots=="2d"){

var g2d_tip = [];
g2d_tip[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"g2dtooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 100px;height: 50px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;

d3.selectAll("."+div_id+"-"+"g2dplot").data(g2dData, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.selectAll("paths")
.data(function(d,i){
var sc_data = [];
for(var l=0;l<d.bar.length;l++){
var objj = {"br":d.bar[l],"nm":d.label[l],"cte":d.category[l]};
sc_data.push(objj);
}
return sc_data;
})
.enter()
.append("path")
.attr("d",d3.symbol().size(loci_width*loci_width).type(gnbox))
.attr("transform",function (d) { return "translate("+(loci_width/2)+","+g2dy(d.br)+")"; })
.style("fill", function(d){
//if(d.br >= 0){ return plot_color;} else {return "red";}
//var scatter_fill;
switch(plot_f[0]){
case "eq":
if(d.br == plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gt":
if(d.br > plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gte":
if(d.br >= plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "lt":
if(d.br < plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "lte":
if(d.br <= plot_f[1]){return plot_f[2]; } else {return plot_color;}
case "gtalt":
if(d.br > plot_f[1] && d.br < plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gtealte":
if(d.br >= plot_f[1] && d.br <= plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gtolt":
if(d.br > plot_f[1] || d.br < plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "gteolte":
if(d.br >= plot_f[1] || d.br <= plot_f[2]){return plot_f[3]; } else {return plot_color;}
case "none":
return plot_color;
case "col":
if(plot_f[1] == "byCategory"){
return scatter_color_map_fn(d.cte);
} else {
return colors[times](d.cte);
}
default:
return plot_color;

}

})
.on("mouseover", function(d) {
g2d_tip[ploidy-1].transition().duration(200).style("opacity", .9).style("background-color","#F5F5F5").style("width",100+(d.nm.length*3)+'px');
g2d_tip[ploidy-1].html("name: "+d.nm+"<br> value: "+d.cte+"<br> position: "+d.br).style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px");
}).on("mouseout", function(d) {	g2d_tip[ploidy-1].transition() .delay(1000).duration(500)	.style("opacity", 0);});
if(ref_line){
for(i=0;i< nLoci.length;i++){

svg.select("#"+div_id+"-"+"g2d-axis-"+nLoci[i].name+"-"+ploidy).append("line")
.attr("x1", 0)
.attr("x2", left_margin +  (win_scale*loci_width))
.attr("y1",  g2dy(refl_pos))
.attr("y2", g2dy(refl_pos))
.attr("stroke", refl_color)
.attr("stroke-dasharray","5,5")
.attr("stroke-width",refl_stroke_w)
.attr("class",div_id+"-"+"overlines2");


}}
}
/**/

if(plots=="tags"){
d3.selectAll("."+div_id+"-"+"tags").data(tagData, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.append("line")
.attr("x1", loci_width/2)
.attr("x2", loci_width/2)
.attr("y1", 0 )
.attr("y2", plot_height )
.attr("stroke", "black");

d3.selectAll("."+div_id+"-"+"tags").data(tagData, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.append("circle")
.attr("cx", loci_width/2)
.attr("cy", loci_width/2)
.attr("r", loci_width/2 )
.attr("fill", tagColor);

}

if(plots=="line"){

  // global domain of data 

  global_domain_line = [];

  for(var i = 0;i < chDataReducedLineLoci.length;i++){

    line_data = chDataReducedLineLoci[i].value; 

    global_domain_line.push(d3.extent(line_data, function(d) { return d.value} ))

  }

  print(global_domain_line)

  if(parseFloat(plot_y_domain[0]) === 0 && parseFloat(plot_y_domain[1]) === 0){
    var line_y = d3.scaleLinear()
    .domain(d3.extent(global_domain_line.flat()))
    .range([ plot_height,0]).nice();
    } else {
    var line_y = d3.scaleLinear()
    .domain(plot_y_domain)
    .range([ plot_height,0]).nice();
    }

  for(var i = 0;i < chDataReducedLineLoci.length;i++){

    line_data = chDataReducedLineLoci[i].value;

    

    d3.select("#line-"+ chDataReducedLineLoci[i].key)
    .append("path")
      .datum(line_data.sort((a, b) => b.win-a.win))
      .attr("fill", "none")
      .attr("stroke", "steelblue")
      .attr("stroke-width", 1.5)
      .attr("d", d3.line()
        .x(function(d) { return (d.win-1)*loci_width+(loci_width/2) })
        .y(function(d) { return line_y(d.value) })
        )


  }



}




} else {
if(heat_scale=="categorical"){
var chDataReduced = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a=[];for(var i=0; i< v.length;i++){ a.push(v[i].data);   };return a[0];})
.entries(chData);


var chDataReducedBarData = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].data);};return bar;})
.entries(chData);
var chDataReducedBarLabel = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].name);};return bar;})
.entries(chData);

var chDataReducedBarLink = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { v=v.sort(function(a,b){return Math.abs(a['ch_start']) - Math.abs(b['ch_start'])});bar=[]; for(var g=0;g<v.length;g++){bar.push(v[g].hlink);};return bar;})
.entries(chData);
var chDataReducedCount = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) { return v.length;})
.entries(chData);
var chDataRange = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a,b;for(var i=0; i< v.length;i++){ a = v[i].loci_start; b=v[i].loci_end;    };return ""+a+"-"+b+"bp ";})
.entries(chData);

for(var i = 0;i < chDataReduced.length;i++) {


chDataReduced[i].range = chDataRange[i].value;

chDataReduced[i].bar = chDataReducedBarData[i].value;
chDataReduced[i].label = chDataReducedBarLabel[i].value;
chDataReduced[i].count = chDataReducedCount[i].value;
chDataReduced[i].hlink = chDataReducedBarLink[i].value;
}

delete chDataReducedBarData;
delete chDataReducedBarLabel;
delete chDataReducedCount;
delete chDataReducedBarLabel;

var tip = [];
tip[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"chtooltip")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 180px;height: 140px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;

var tip2 = [];
tip2[ploidy-1] = d3.select("body").append("div")
.attr("class", div_id+"-"+"chtooltip2")
.style("opacity", 0)
.attr("style", "position: absolute;text-align: center;width: 180px;height: 140px;	 padding: 2px;	font: 12px sans-serif;	box-sizing: border-box;	border: 0px;border-radius: 8px;pointer-events: none;") ;
var clickFlag=false;
if(interactivity){

d3.selectAll("."+div_id+"-"+"chLoc").data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");}).on("mouseover", function(d) {tip[ploidy-1].transition().duration(200)	.style("opacity", .9).style("background-color","#F5F5F5").style("width",180+d.count*(d3.max(d.label).length*3)+'px');
tip[ploidy-1].html("<div  style='float:top;position:relative;height:35%;' ><div style='border-radius: 5px;float:left;position:relative;height:55px;width:30%;background-color:"+colors[times](d.value)+
";' ></div><div  style='float:left;position:relative;width:70%;' ><p ><font size='1' color='grey' >  </font> <font size='1'  ><b> "+""+
" <br><font size='1' color='grey' > Count:</font>"+d.count+" <br><font size='1' color='grey' >  </font> "+""+
"</font><br><font size='1' color='grey' >  </font><font size='1'>"+""+"</p></font></div></div><div id ='"+div_id+"-microBar"+ploidy+"' style='float:top;position:relative;height:35%'></div><br><div  style='float:top;height:20%;position:relative;color:'black"+
";'><font size='2' >"+d.range+"</font></p></div>")
.style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px"); var svgWidth = 175+d.count*(d3.max(d.label).length)*3;var svgHeight = 50;var barPadding = 5;var barWidth = (svgWidth / d.bar.length);var mysvg = d3.select('#'+div_id+'-microBar'+ploidy).append("svg").attr("width",svgWidth).attr("heigth",svgHeight);var barChart = mysvg.selectAll('rect').data(d.bar).enter().append('rect').attr('height',svgHeight-40).attr('width', barWidth).attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';}).style("fill",function(d){return colors[times](d);});
t_data=[];
for(j=0;j<d.label.length;j++){
t_data.push({lb:d.label[j],hl:d.hlink[j]});
}

var ttext=mysvg.selectAll(".node").data(t_data).enter().append("svg:a").attr("xlink:href", function(d,i){ return d.hl;})
.append("svg:text").attr("y",svgHeight*0.4)
.attr("transform","rotate(90)")
.attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';})
.attr("font-size", "9px").style("pointer-events", "all")
.attr("class","ttt").style("cursor","pointer")
.text(function(d,i){return d.lb;});  })
.on("mouseout", function(d) {	tip[ploidy-1].transition() .delay(1000).duration(500)	.style("opacity", 0);}).style("visibility","visible").style("fill",function(d){ return colors[times](d.value);})
.on("click", function(d) {


if(clickFlag){
tip2[ploidy-1].style("opacity", 0);
}else{
tip2[ploidy-1].transition().duration(200)	.style("opacity", .9).style("background-color","#F5F5F5").style("width",180+d.count*(d3.max(d.label).length*3)+'px');
tip2[ploidy-1].html("<div  style='float:top;position:relative;height:35%;' ><div style='border-radius: 5px;float:left;position:relative;height:55px;width:30%;background-color:"+colors[times](d.value)+
";' ></div><div  style='float:left;position:relative;width:70%;' ><p ><font size='1' color='grey' >  </font> <font size='1'  ><b> "+""+
" <br><font size='1' color='grey' > Count:</font>"+d.count+" <br><font size='1' color='grey' >  </font> "+""+
"</font><br><font size='1' color='grey' >  </font><font size='1'>"+""+"</p></font></div></div><div id ='"+div_id+"-microBar2"+ploidy+"' style='float:top;position:relative;height:35%'></div><br><div  style='float:top;height:20%;position:relative;color:'black"+
";'><font size='2' >"+d.range+"</font></p></div>")
.style("left", (d3.event.pageX) + "px")
.style("top", (d3.event.pageY - 18) + "px"); var svgWidth = 175+d.count*(d3.max(d.label).length)*3;var svgHeight = 50;var barPadding = 5;var barWidth = (svgWidth / d.bar.length);var mysvg2 = d3.select('#'+div_id+'-microBar2'+ploidy).append("svg").attr("width",svgWidth).attr("heigth",svgHeight);var barChart = mysvg2.selectAll('rect').data(d.bar).enter().append('rect').attr('height',svgHeight-40).attr('width', barWidth).attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';}).style("fill",function(d){return colors[times](d);});
t_data=[];
for(j=0;j<d.label.length;j++){
t_data.push({lb:d.label[j],hl:d.hlink[j]});
}

var ttext=mysvg2.selectAll(".node").data(t_data).enter().append("svg:a").attr("xlink:href", function(d,i){ return d.hl;})
.append("svg:text").attr("y",svgHeight*0.4)
.attr("transform","rotate(90)")
.attr('transform', function (d, i) { var translate = [barWidth * i, 0];  return 'translate('+ translate +')';})
.attr("font-size", "9px").style("pointer-events", "all")
.attr("class","ttt").style("cursor","pointer")
.text(function(d,i){return d.lb;}); } return clickFlag = !clickFlag; });

} else {
//without interactvity
d3.selectAll("."+div_id+"-"+"chLoc")
.data(chDataReduced, function(d) {return (d && d.key)||d3.select(this).attr("id");})
.style("visibility","visible")
.style("fill",function(d){ return colors[times](d.value);});


//witout interactivity end
}



if(labels){

var chDataLabels = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a=[];for(var i=0; i< v.length;i++){ a.push(v[i].name);   };return a[0];})
.entries(chData);
var chDataLabel = d3.nest()
.key(function(d) { return d.loci; })
.rollup(function(v) {  var a=[];for(var i=0; i< v.length;i++){ a.push(v[i].label);   };return a[0];})
.entries(chData);
for(var i = 0;i < chDataReduced.length;i++) {


chDataReduced[i].label_id = chDataLabel[i].value;
chDataReduced[i].labels =chDataLabels[i].value;

}

d3.selectAll("."+div_id+"-"+"labels").data(chDataReduced, function(d) {return (d && d.label_id)||d3.select(this).attr("id");})
.text(function(d){ return d.labels;})
.attr("transform", function (d) {
var xRot = d3.select(this).attr("x");
var yRot = d3.select(this).attr("y");
return `rotate(${label_angle}, ${xRot},  ${yRot} )`
});
}
}

}













//heatmap code end
}


}




// Global variable for chromoMap 

var chromoMap = new chromoMapPlot();


function saveSvg2(idd, name) {
  svgEl = document.getElementById(idd);
  svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");

  var rect = document.createElementNS("http://www.w3.org/2000/svg", 'rect');
rect.setAttribute('x', '0');
rect.setAttribute('y', '0');
rect.setAttribute('height', '100%');
rect.setAttribute('width', '100%');
rect.setAttribute('fill', "white");
svgEl.prepend(rect);

  var svgData = svgEl.outerHTML;
  var preface = '<?xml version="1.0" standalone="no"?>\r\n';
  var svgBlob = new Blob([preface, svgData], {type:"image/svg+xml;charset=utf-8"});
  var svgUrl = URL.createObjectURL(svgBlob);
  var downloadLink = document.createElement("a");
  downloadLink.href = svgUrl;
  downloadLink.download = name;
  document.body.appendChild(downloadLink);
  downloadLink.click();
  document.body.removeChild(downloadLink);


}


(function() {
  const out$ = typeof exports != 'undefined' && exports || typeof define != 'undefined' && {} || this || window;
  if (typeof define !== 'undefined') define('save-svg-as-png', [], () => out$);
  out$.default = out$;

  const xmlNs = 'http://www.w3.org/2000/xmlns/';
  const xhtmlNs = 'http://www.w3.org/1999/xhtml';
  const svgNs = 'http://www.w3.org/2000/svg';
  const doctype = '<?xml version="1.0" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd" [<!ENTITY nbsp "&#160;">]>';
  const urlRegex = /url\(["']?(.+?)["']?\)/;
  const fontFormats = {
    woff2: 'font/woff2',
    woff: 'font/woff',
    otf: 'application/x-font-opentype',
    ttf: 'application/x-font-ttf',
    eot: 'application/vnd.ms-fontobject',
    sfnt: 'application/font-sfnt',
    svg: 'image/svg+xml'
  };

  const isElement = obj => obj instanceof HTMLElement || obj instanceof SVGElement;
  const requireDomNode = el => {
    if (!isElement(el)) throw new Error(`an HTMLElement or SVGElement is required; got ${el}`);
  };
  const requireDomNodePromise = el =>
    new Promise((resolve, reject) => {
      if (isElement(el)) resolve(el)
      else reject(new Error(`an HTMLElement or SVGElement is required; got ${el}`));
    })
  const isExternal = url => url && url.lastIndexOf('http',0) === 0 && url.lastIndexOf(window.location.host) === -1;

  const getFontMimeTypeFromUrl = fontUrl => {
    const formats = Object.keys(fontFormats)
      .filter(extension => fontUrl.indexOf(`.${extension}`) > 0)
      .map(extension => fontFormats[extension]);
    if (formats) return formats[0];
    console.error(`Unknown font format for ${fontUrl}. Fonts may not be working correctly.`);
    return 'application/octet-stream';
  };

  const arrayBufferToBase64 = buffer => {
    let binary = '';
    const bytes = new Uint8Array(buffer);
    for (let i = 0; i < bytes.byteLength; i++) binary += String.fromCharCode(bytes[i]);
    return window.btoa(binary);
  }

  const getDimension = (el, clone, dim) => {
    const v =
      (el.viewBox && el.viewBox.baseVal && el.viewBox.baseVal[dim]) ||
      (clone.getAttribute(dim) !== null && !clone.getAttribute(dim).match(/%$/) && parseInt(clone.getAttribute(dim))) ||
      el.getBoundingClientRect()[dim] ||
      parseInt(clone.style[dim]) ||
      parseInt(window.getComputedStyle(el).getPropertyValue(dim));
    return typeof v === 'undefined' || v === null || isNaN(parseFloat(v)) ? 0 : v;
  };

  const getDimensions = (el, clone, width, height) => {
    if (el.tagName === 'svg') return {
      width: width || getDimension(el, clone, 'width'),
      height: height || getDimension(el, clone, 'height')
    };
    else if (el.getBBox) {
      const {x, y, width, height} = el.getBBox();
      return {
        width: x + width,
        height: y + height
      };
    }
  };

  const reEncode = data =>
    decodeURIComponent(
      encodeURIComponent(data)
        .replace(/%([0-9A-F]{2})/g, (match, p1) => {
          const c = String.fromCharCode(`0x${p1}`);
          return c === '%' ? '%25' : c;
        })
    );

  const uriToBlob = uri => {
    const byteString = window.atob(uri.split(',')[1]);
    const mimeString = uri.split(',')[0].split(':')[1].split(';')[0]
    const buffer = new ArrayBuffer(byteString.length);
    const intArray = new Uint8Array(buffer);
    for (let i = 0; i < byteString.length; i++) {
      intArray[i] = byteString.charCodeAt(i);
    }
    return new Blob([buffer], {type: mimeString});
  };

  const query = (el, selector) => {
    if (!selector) return;
    try {
      return el.querySelector(selector) || el.parentNode && el.parentNode.querySelector(selector);
    } catch(err) {
      console.warn(`Invalid CSS selector "${selector}"`, err);
    }
  };

  const detectCssFont = (rule, href) => {
    // Match CSS font-face rules to external links.
    // @font-face {
    //   src: local('Abel'), url(https://fonts.gstatic.com/s/abel/v6/UzN-iejR1VoXU2Oc-7LsbvesZW2xOQ-xsNqO47m55DA.woff2);
    // }
    const match = rule.cssText.match(urlRegex);
    const url = (match && match[1]) || '';
    if (!url || url.match(/^data:/) || url === 'about:blank') return;
    const fullUrl =
      url.startsWith('../') ? `${href}/../${url}`
      : url.startsWith('./') ? `${href}/.${url}`
      : url;
    return {
      text: rule.cssText,
      format: getFontMimeTypeFromUrl(fullUrl),
      url: fullUrl
    };
  };

  const inlineImages = el => Promise.all(
    Array.from(el.querySelectorAll('image')).map(image => {
      let href = image.getAttributeNS('http://www.w3.org/1999/xlink', 'href') || image.getAttribute('href');
      if (!href) return Promise.resolve(null);
      if (isExternal(href)) {
        href += (href.indexOf('?') === -1 ? '?' : '&') + 't=' + new Date().valueOf();
      }
      return new Promise((resolve, reject) => {
        const canvas = document.createElement('canvas');
        const img = new Image();
        img.crossOrigin = 'anonymous';
        img.src = href;
        img.onerror = () => reject(new Error(`Could not load ${href}`));
        img.onload = () => {
          canvas.width = img.width;
          canvas.height = img.height;
          canvas.getContext('2d').drawImage(img, 0, 0);
          image.setAttributeNS('http://www.w3.org/1999/xlink', 'href', canvas.toDataURL('image/png'));
          resolve(true);
        };
      });
    })
  );

  const cachedFonts = {};
  const inlineFonts = fonts => Promise.all(
    fonts.map(font =>
      new Promise((resolve, reject) => {
        if (cachedFonts[font.url]) return resolve(cachedFonts[font.url]);

        const req = new XMLHttpRequest();
        req.addEventListener('load', () => {
          // TODO: it may also be worth it to wait until fonts are fully loaded before
          // attempting to rasterize them. (e.g. use https://developer.mozilla.org/en-US/docs/Web/API/FontFaceSet)
          const fontInBase64 = arrayBufferToBase64(req.response);
          const fontUri = font.text.replace(urlRegex, `url("data:${font.format};base64,${fontInBase64}")`)+'\n';
          cachedFonts[font.url] = fontUri;
          resolve(fontUri);
        });
        req.addEventListener('error', e => {
          console.warn(`Failed to load font from: ${font.url}`, e);
          cachedFonts[font.url] = null;
          resolve(null);
        });
        req.addEventListener('abort', e => {
          console.warn(`Aborted loading font from: ${font.url}`, e);
          resolve(null);
        });
        req.open('GET', font.url);
        req.responseType = 'arraybuffer';
        req.send();
      })
    )
  ).then(fontCss => fontCss.filter(x => x).join(''));

  let cachedRules = null;
  const styleSheetRules = () => {
    if (cachedRules) return cachedRules;
    return cachedRules = Array.from(document.styleSheets).map(sheet => {
      try {
        return {rules: sheet.cssRules, href: sheet.href};
      } catch (e) {
        console.warn(`Stylesheet could not be loaded: ${sheet.href}`, e);
        return {};
      }
    });
  };

  const inlineCss = (el, options) => {
    const {
      selectorRemap,
      modifyStyle,
      modifyCss,
      fonts,
      excludeUnusedCss
    } = options || {};
    const generateCss = modifyCss || ((selector, properties) => {
      const sel = selectorRemap ? selectorRemap(selector) : selector;
      const props = modifyStyle ? modifyStyle(properties) : properties;
      return `${sel}{${props}}\n`;
    });
    const css = [];
    const detectFonts = typeof fonts === 'undefined';
    const fontList = fonts || [];
    styleSheetRules().forEach(({rules, href}) => {
      if (!rules) return;
      Array.from(rules).forEach(rule => {
        if (typeof rule.style != 'undefined') {
          if (query(el, rule.selectorText)) css.push(generateCss(rule.selectorText, rule.style.cssText));
          else if (detectFonts && rule.cssText.match(/^@font-face/)) {
            const font = detectCssFont(rule, href);
            if (font) fontList.push(font);
          } else if (!excludeUnusedCss) {
            css.push(rule.cssText);
          }
        }
      });
    });

    return inlineFonts(fontList).then(fontCss => css.join('\n') + fontCss);
  };

  const downloadOptions = () => {
    if (!navigator.msSaveOrOpenBlob && !('download' in document.createElement('a'))) {
      return {popup: window.open()};
    }
  };

  out$.prepareSvg = (el, options, done) => {
    requireDomNode(el);
    const {
      left = 0,
      top = 0,
      width: w,
      height: h,
      scale = 1,
      responsive = false,
      excludeCss = false,
    } = options || {};

    return inlineImages(el).then(() => {
      let clone = el.cloneNode(true);
      clone.style.backgroundColor = (options || {}).backgroundColor || el.style.backgroundColor;
      const {width, height} = getDimensions(el, clone, w, h);

      if (el.tagName !== 'svg') {
        if (el.getBBox) {
          if (clone.getAttribute('transform') != null) {
            clone.setAttribute('transform', clone.getAttribute('transform').replace(/translate\(.*?\)/, ''));
          }
          const svg = document.createElementNS('http://www.w3.org/2000/svg','svg');
          svg.appendChild(clone);
          clone = svg;
        } else {
          console.error('Attempted to render non-SVG element', el);
          return;
        }
      }

      clone.setAttribute('version', '1.1');
      clone.setAttribute('viewBox', [left, top, width, height].join(' '));
      if (!clone.getAttribute('xmlns')) clone.setAttributeNS(xmlNs, 'xmlns', svgNs);
      if (!clone.getAttribute('xmlns:xlink')) clone.setAttributeNS(xmlNs, 'xmlns:xlink', 'http://www.w3.org/1999/xlink');

      if (responsive) {
        clone.removeAttribute('width');
        clone.removeAttribute('height');
        clone.setAttribute('preserveAspectRatio', 'xMinYMin meet');
      } else {
        clone.setAttribute('width', width * scale);
        clone.setAttribute('height', height * scale);
      }

      Array.from(clone.querySelectorAll('foreignObject > *')).forEach(foreignObject => {
        foreignObject.setAttributeNS(xmlNs, 'xmlns', foreignObject.tagName === 'svg' ? svgNs : xhtmlNs);
      });

      if (excludeCss) {
        const outer = document.createElement('div');
        outer.appendChild(clone);
        const src = outer.innerHTML;
        if (typeof done === 'function') done(src, width, height);
        else return {src, width, height};
      } else {
        return inlineCss(el, options).then(css => {
          const style = document.createElement('style');
          style.setAttribute('type', 'text/css');
          style.innerHTML = `<![CDATA[\n${css}\n]]>`;

          const defs = document.createElement('defs');
          defs.appendChild(style);
          clone.insertBefore(defs, clone.firstChild);

          const outer = document.createElement('div');
          outer.appendChild(clone);
          const src = outer.innerHTML.replace(/NS\d+:href/gi, 'xmlns:xlink="http://www.w3.org/1999/xlink" xlink:href');

          if (typeof done === 'function') done(src, width, height);
          else return {src, width, height};
        });
      }
    });
  };

  out$.svgAsDataUri = (el, options, done) => {
    requireDomNode(el);
    return out$.prepareSvg(el, options)
      .then(({src, width, height}) => {
          const svgXml = `data:image/svg+xml;base64,${window.btoa(reEncode(doctype+src))}`;
          if (typeof done === 'function') {
              done(svgXml, width, height);
          }
          return svgXml;
      });
  };

  out$.svgAsPngUri = (el, options, done) => {
    requireDomNode(el);
    const {
      encoderType = 'image/png',
      encoderOptions = 0.8,
      canvg
    } = options || {};

    const convertToPng = ({src, width, height}) => {
      const canvas = document.createElement('canvas');
      const context = canvas.getContext('2d');
      const pixelRatio = window.devicePixelRatio || 1;

      canvas.width = width * pixelRatio;
      canvas.height = height * pixelRatio;
      canvas.style.width = `${canvas.width}px`;
      canvas.style.height = `${canvas.height}px`;
      context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);

      if (canvg) canvg(canvas, src);
      else context.drawImage(src, 0, 0);

      let png;
      try {
        png = canvas.toDataURL(encoderType, encoderOptions);
      } catch (e) {
        if ((typeof SecurityError !== 'undefined' && e instanceof SecurityError) || e.name === 'SecurityError') {
          console.error('Rendered SVG images cannot be downloaded in this browser.');
          return;
        } else throw e;
      }
      if (typeof done === 'function') done(png, canvas.width, canvas.height);
      return Promise.resolve(png);
    }

    if (canvg) return out$.prepareSvg(el, options).then(convertToPng);
    else return out$.svgAsDataUri(el, options).then(uri => {
      return new Promise((resolve, reject) => {
        const image = new Image();
        image.onload = () => resolve(convertToPng({
          src: image,
          width: image.width,
          height: image.height
        }));
        image.onerror = () => {
          reject(`There was an error loading the data URI as an image on the following SVG\n${window.atob(uri.slice(26))}Open the following link to see browser's diagnosis\n${uri}`);
        }
        image.src = uri;
      })
    });
  };

  out$.download = (name, uri, options) => {
    if (navigator.msSaveOrOpenBlob) navigator.msSaveOrOpenBlob(uriToBlob(uri), name);
    else {
      const saveLink = document.createElement('a');
      if ('download' in saveLink) {
        saveLink.download = name;
        saveLink.style.display = 'none';
        document.body.appendChild(saveLink);
        try {
          const blob = uriToBlob(uri);
          const url = URL.createObjectURL(blob);
          saveLink.href = url;
          saveLink.onclick = () => requestAnimationFrame(() => URL.revokeObjectURL(url));
        } catch (e) {
          console.error(e);
          console.warn('Error while getting object URL. Falling back to string URL.');
          saveLink.href = uri;
        }
        saveLink.click();
        document.body.removeChild(saveLink);
      } else if (options && options.popup) {
        options.popup.document.title = name;
        options.popup.location.replace(uri);
      }
    }
  };

  out$.saveSvg = (el, name, options) => {
    const downloadOpts = downloadOptions(); // don't inline, can't be async
    return requireDomNodePromise(el)
      .then(el => out$.svgAsDataUri(el, options || {}))
      .then(uri => out$.download(name, uri, downloadOpts));
  };

  out$.saveSvgAsPng = (el, name, options) => {
    const downloadOpts = downloadOptions(); // don't inline, can't be async
    return requireDomNodePromise(el)
      .then(el => out$.svgAsPngUri(el, options || {}))
      .then(uri => out$.download(name, uri, downloadOpts));
  };
})();

