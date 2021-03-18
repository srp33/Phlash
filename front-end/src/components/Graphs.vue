<template>
  <div id="columns">
    <div class="column">
      <h4>Direct Sequences</h4>
      <Plotly
        class="graph"
        v-if="frame == 1"
        :data="data1"
        :layout="layout1Selected"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame !== 1"
        :data="data1"
        :layout="layout1"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame == 2"
        :data="data2"
        :layout="layout2Selected"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame !== 2"
        :data="data2"
        :layout="layout2"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame === 3"
        :data="data3"
        :layout="layout3Selected"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame !== 3"
        :data="data3"
        :layout="layout3"
        :display-mode-bar="false"
      ></Plotly>
    </div>
    <div class="column">
      <h4>Complementary Sequences</h4>
      <Plotly
        class="graph"
        v-if="frame === 4"
        :data="data4"
        :layout="layout4Selected"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame !== 4"
        :data="data4"
        :layout="layout4"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame === 5"
        :data="data5"
        :layout="layout5Selected"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame !== 5"
        :data="data5"
        :layout="layout5"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame === 6"
        :data="data6"
        :layout="layout6Selected"
        :display-mode-bar="false"
      ></Plotly>
      <Plotly
        class="graph"
        v-if="frame !== 6"
        :data="data6"
        :layout="layout6"
        :display-mode-bar="false"
      ></Plotly>
    </div>
  </div>
</template>

<script>
import { Plotly } from 'vue-plotly';
export default {
  name: 'Graphs',
  components: {
    Plotly,
  },
  props: {
    frame: Number,
    left: Number,
    right: Number,
    data1: Array,
    data2: Array,
    data3: Array,
    data4: Array,
    data5: Array,
    data6: Array,
    prevRight: Number,
    nextLeft: Number,
  },

  data() {
    return {
      // ----- LAYOUT 1 -----
      layout1: {
        autosize: false,
        width: 590,
        height: 200,
        color: 'black',
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          dtick: 200,
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 1',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
        ],
        annotations: [
          {
            x: this.left,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>left',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
          {
            x: this.right,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>right',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
        ],
      },
      // ----- LAYOUT 2 -----
      layout2: {
        autosize: false,
        width: 590,
        height: 200,
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          dtick: 200,
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 2',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
        ],
      },
      // ----- LAYOUT 3 -----
      layout3: {
        autosize: false,
        width: 590,
        height: 225,
        margin: {
          l: 52,
          r: 50,
          b: 60,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          title: 'Nucleotide Position',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 3',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
        ],
      },
      // ----- LAYOUT 4 -----
      layout4: {
        autosize: false,
        width: 590,
        height: 200,
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 4',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
        ],
        annotations: [
          {
            x: this.left,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>left',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
          {
            x: this.right,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>right',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
        ],
      },
      // ----- LAYOUT 5 -----
      layout5: {
        autosize: false,
        width: 590,
        height: 200,
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 5',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
        ],
      },
      // ----- LAYOUT 6 -----
      layout6: {
        autosize: false,
        width: 590,
        height: 225,
        margin: {
          l: 52,
          r: 50,
          b: 60,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          title: 'Nucleotide Position',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 6',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
        ],
      },
      // ----- LAYOUT 1 Selected -----
      layout1Selected: {
        autosize: false,
        paper_bgcolor: 'rgba(219,219,219,0.5)',
        width: 590,
        height: 200,
        color: 'black',
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },

        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          dtick: 200,
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 1',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.left,
            y0: 0.75,
            x1: this.right,
            yref: 'paper',
            y1: 0.75,
            line: {
              color: '#1b9e77',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.prevRight,
            y0: 0,
            x1: this.prevRight,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.nextLeft,
            y0: 0,
            x1: this.nextLeft,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
        ],
        annotations: [
          {
            x: this.left,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>left',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
          {
            x: this.right,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>right',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
        ],
      },
      // ----- LAYOUT 2 Selected -----
      layout2Selected: {
        autosize: false,
        paper_bgcolor: 'rgba(219,219,219,0.5)',
        width: 590,
        height: 200,
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          dtick: 200,
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 2',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.left,
            y0: 0.75,
            x1: this.right,
            yref: 'paper',
            y1: 0.75,
            line: {
              color: '#1b9e77',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.prevRight,
            y0: 0,
            x1: this.prevRight,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.nextLeft,
            y0: 0,
            x1: this.nextLeft,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
        ],
      },
      // ----- LAYOUT 3 Selected -----
      layout3Selected: {
        autosize: false,
        paper_bgcolor: 'rgba(219,219,219,0.5)',
        width: 590,
        height: 225,
        margin: {
          l: 52,
          r: 50,
          b: 60,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          title: 'Nucleotide Position',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 3',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.left,
            y0: 0.75,
            x1: this.right,
            yref: 'paper',
            y1: 0.75,
            line: {
              color: '#1b9e77',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.prevRight,
            y0: 0,
            x1: this.prevRight,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.nextLeft,
            y0: 0,
            x1: this.nextLeft,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
        ],
      },
      // ----- LAYOUT 4 Selected -----
      layout4Selected: {
        autosize: false,
        paper_bgcolor: 'rgba(219,219,219,0.5)',
        width: 590,
        height: 200,
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 4',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.left,
            y0: 0.75,
            x1: this.right,
            yref: 'paper',
            y1: 0.75,
            line: {
              color: '#1b9e77',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.prevRight,
            y0: 0,
            x1: this.prevRight,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.nextLeft,
            y0: 0,
            x1: this.nextLeft,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
        ],
        annotations: [
          {
            x: this.left,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>left',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
          {
            x: this.right,
            y: 1,
            xref: 'x',
            yref: 'y',
            text: 'current<br>right',
            font: {
              family: 'Avenir, Helvetica, Arial, sans-serif',
              size: 11,
              color: 'black',
            },
            showarrow: true,
            ax: 0,
            ay: -20,
          },
        ],
      },
      // ----- LAYOUT 5 Selected -----
      layout5Selected: {
        autosize: false,
        paper_bgcolor: 'rgba(219,219,219,0.5)',
        width: 590,
        height: 200,
        margin: {
          l: 52,
          r: 50,
          b: 30,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 5',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.left,
            y0: 0.75,
            x1: this.right,
            yref: 'paper',
            y1: 0.75,
            line: {
              color: '#1b9e77',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.prevRight,
            y0: 0,
            x1: this.prevRight,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.nextLeft,
            y0: 0,
            x1: this.nextLeft,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
        ],
      },
      // ----- LAYOUT 6 Selected -----
      layout6Selected: {
        autosize: false,
        paper_bgcolor: 'rgba(219,219,219,0.5)',
        width: 590,
        height: 225,
        margin: {
          l: 52,
          r: 50,
          b: 60,
          t: 30,
          pad: 4,
        },
        spikedistance: 200,
        hoverdistance: 100,
        xaxis: {
          title: 'Nucleotide Position',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          dtick: 200,
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
          showcrossline: true,
          spikemode: 'across',
          spikecolor: '#1b9e77',
          spikedash: 'solid',
          spikethickness: 2,
        },
        yaxis: {
          title: 'Frame 6',
          titlefont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 14,
            color: 'grey',
          },
          showgrid: false,
          showline: true,
          linewidth: 2,
          ticks: 'outside',
          ticklen: 7,
          tickwidth: 2,
          range: [0, 1],
          dtick: 0.5,
          tickformat: '.1f',
          tickfont: {
            family: 'Avenir, Helvetica, Arial, sans-serif',
            size: 12,
            color: 'black',
          },
        },
        shapes: [
          {
            type: 'line',
            x0: this.left,
            y0: 0,
            x1: this.left,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.right,
            y0: 0,
            x1: this.right,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#d95f02',
              width: 1.5,
              dash: 'dot',
            },
          },
          {
            type: 'line',
            x0: this.left,
            y0: 0.75,
            x1: this.right,
            yref: 'paper',
            y1: 0.75,
            line: {
              color: '#1b9e77',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.prevRight,
            y0: 0,
            x1: this.prevRight,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
          {
            type: 'line',
            x0: this.nextLeft,
            y0: 0,
            x1: this.nextLeft,
            yref: 'paper',
            y1: 1,
            line: {
              color: '#7570b3',
              width: 1.0,
            },
          },
        ],
      },
    };
  },
};
</script>

<style scoped>
.columns:after {
  display: table;
  clear: both;
}

.column {
  width: 50%;
  float: left;
}

/* --------------------- */
/* Responsive Design */
@media only screen and (max-width: 1200px) {
  .column {
    width: 100%;
    margin-bottom: 30px;
  }
}
</style>