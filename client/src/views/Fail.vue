<template>
<div class="container">
   <div class="headers">
      <h1>ID: {{ $route.params.id }}</h1>
      <h4>Start: {{ currentCDS.start }}</h4>
      <h4>Stop: {{ currentCDS.stop }}</h4>
      <h4>Status: {{ currentCDS.status }}</h4>
   </div>
   <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>
         Choose a new start for this gene call based on the information given
         below.<br/>Here are your other possible start positions: 
         <span v-for="(curr, index) in startOptions" :key=index>
            <span v-if="index===startOptions.length-1">{{ curr.start_position }}.</span>
            <span v-else>{{ curr.start_position }}, </span>
         </span>
      </p>
      <p><strong>Your selected start position:</strong> {{ newStart }}</p>
      <p><strong>Your selected function:</strong> {{ newFunction }}</p>
      <button type="button" class="btn btn-light" @click="editCDS" v-if="newFunction !== 'None' && newStart !== 0">
         <strong>Update {{ $route.params.id }}</strong>
      </button>
   </div>
   <div class="columns-wrapper">
      <h3 style="text-align: center; margin: 20px;">Coding Potential</h3>
      <div class="subheader">
         <!-- <p>
            The graph below displays the coding potential in each open
            reading frame (ORF). The first three plots represent direct
            sequences, and the latter three represent complementary
            sequences. See the table below the graph for more information.
         </p>
         <p>
            <strong>Instructions:</strong> Hover your mouse over the graph to
            see the exact base positions. Then, scroll down to the table
            below and <em>select a new start position, if
            appropriate</em>.
         </p> -->
      </div>
   </div>
   <div class="coding-potential-table">
      <div class="subheader">
         <p>
            This table displays the average coding potential per
            frame between the specified start position and {{
            currentCDS.stop }} (the graph above is a visual
            representation of this data). Frames are counted 1-6
            (direct 1-3 and complementary 4-6). Select a new start
            position, if appropriate.
         </p>
      </div>
      <div class="table-responsive">
         <table id="cp-table" class="table table-hover">
            <thead>
               <tr>
                  <th scope="col">Start Position</th>
                  <th scope="col">Frame 1</th>
                  <th scope="col">Frame 2</th>
                  <th scope="col">Frame 3</th>
                  <th scope="col">Frame 4</th>
                  <th scope="col">Frame 5</th>
                  <th scope="col">Frame 6</th>
                  <th scope="col">Action</th>
               </tr>
            </thead>
            <tbody>
               <tr v-for="(curr, index) in startOptions" :key="index">
                  <th>{{ curr.start_position }}</th>
                  <td>{{ curr.frame_1 }}</td>
                  <td>{{ curr.frame_2 }}</td>
                  <td>{{ curr.frame_3 }}</td>
                  <td>{{ curr.frame_4 }}</td>
                  <td>{{ curr.frame_5 }}</td>
                  <td>{{ curr.frame_6 }}</td>
                  <td>
                     <button class="btn btn-dark btn-sm" @click="setStart(curr.start_position)">
                        <strong>Select</strong>
                     </button>
                  </td>
               </tr>
            </tbody>
         </table>
      </div>
   </div>
   <!-- <div class="columns">
      <div class="column direct">
         <div id="chart-wrapper">
            <div id="chart-line">
               <apexchart type="line" height="170" :options="options1" :series="series1"></apexchart>
            </div>
            <div id="chart-line2">
               <apexchart type="line" height="170" :options="options2" :series="series2"></apexchart>
            </div>
            <div id="chart-area">
               <apexchart type="line" height="193" :options="options3" :series="series3"></apexchart>
            </div>
         </div>
      </div>
      <div class="column complementary">
         <div id="chart-wrapper">
            <div id="chart-line">
               <apexchart type="line" height="170" :options="options4" :series="series4"></apexchart>
            </div>
            <div id="chart-line2">
               <apexchart type="line" height="170" :options="options5" :series="series5"></apexchart>
            </div>
            <div id="chart-area">
               <apexchart type="line" height="193" :options="options6" :series="series6"></apexchart>
            </div>
         </div>
      </div>
   </div> -->
   <hr/>
   <div class="blast-results">
      <h3 style="text-align: center; margin: 20px;">BLAST Results</h3>
      <div id="accordion">
         <div class="card" v-for="(item, key) in blastResults" :key="key">
            <div class="card-header">
               <h4 class="mb-0">
                  <button class="btn btn-light btn-blast" data-toggle="collapse" aria-expanded="false" 
                     v-bind:data-target="'#'+key" v-bind:aria-controls="key">
                     <strong>{{ key }}</strong>
                  </button>
               </h4>
            </div>
            <div v-bind:id="key" class="collapse" data-parent="#accordion">
               <div class="card-body">
                  <BlastResults :blastResults="item" @newFunction="setFunction" />
               </div>
            </div>
         </div>
      </div>
   </div>
   <div class="info-bottom">
      <p><strong>Your selected start position:</strong> {{ newStart }}</p>
      <p><strong>Your function selection:</strong> {{ newFunction }}</p>
      <button type="button" class="btn btn-light" @click="editCDS" v-if="newFunction !== 'None' && newStart !== 0">
         <strong>Update {{ $route.params.id }}</strong>
      </button>
   </div>
</div>
</template>

<script>
import BlastResults from '../components/BlastResults.vue'
import axios from 'axios'
import VueApexCharts from 'vue-apexcharts'

export default {
   name: 'Fail',
   components: {
      BlastResults
   },
   data () {
      return {
         startOptions: [],
         currentCDS: {
            id: '',
            start: '',
            stop: '',
            strand: '',
            function: '',
            status: '',
         },
         updatedCDS: {
            id: '',
            start: '',
            stop: '',
            strand: '',
            function: '',
            status: '',
         },
         blastResults: [],
         showFunction: false,
         showStart: false,
         newFunction: 'None',
         newStart: 0,
         // series1: [{ data: [] }],
         // series2: [{ data: [] }],
         // series3: [{ data: [] }],
         // series4: [{ data: [] }],
         // series5: [{ data: [] }],
         // series6: [{ data: [] }],
         // options1: {
         //    chart: { id: '1', group: 'direct', type: 'line', height: 170 },
         //    stroke: { width: 1.5 },
         //    colors: ['#008FFB'],
         //    noData: { text: 'Loading...' },
         //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
         //    xaxis: { min: 1, max: 10, tickAmount: 10 }
         // },
         // options2: {
         //    chart: { id: '2', group: 'direct', type: 'line', height: 170 },
         //    stroke: { width: 1.5 },
         //    colors: ['#008FFB'],
         //    noData: { text: 'Loading...' },
         //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 }, 
         //             title: { text: 'Direct Sequences'} },
         //    xaxis: { min: 1, max: 10, tickAmount: 10 }
         // },
         // options3: {
         //    chart: { id: '3', group: 'direct', type: 'line', height: 193 },
         //    stroke: { width: 1.5 },
         //    colors: ['#008FFB'],
         //    noData: { text: 'Loading...' },
         //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
         //    xaxis: { min: 1, max: 10, tickAmount: 10,
         //             title: { text: 'Nucleotide Position'} }
         // },
         // options4: {
         //    chart: { id: '4', group: 'comp', type: 'line', height: 170 },
         //    stroke: { width: 1.5 },
         //    colors: ['#008FFB'],
         //    noData: { text: 'Loading...' },
         //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
         //    xaxis: { min: 1, max: 10, tickAmount: 10 }
         // },
         // options5: {
         //    chart: { id: '5', group: 'comp', type: 'line', height: 170 },
         //    stroke: { width: 1.5 },
         //    colors: ['#008FFB'],
         //    noData: { text: 'Loading...' },
         //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 },
         //             title: { text: 'Complementary Sequences'} },
         //    xaxis: { min: 1, max: 10, tickAmount: 10 }
         // },
         // options6: {
         //    chart: { id: '6', group: 'comp', type: 'line', height: 193 },
         //    stroke: { width: 1.5 },
         //    colors: ['#008FFB'],
         //    noData: { text: 'Loading...' },
         //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
         //    xaxis: { min: 1, max: 10, tickAmount: 10,
         //             title: { text: 'Nucleotide Position'}}
         // },
      }
   },
   created() {
      this.getData(this.$route.params.id);
   },
   methods: {
      getData(cdsID) {
         axios.get(`http://localhost:5000/api/annotations/fail/${this.$route.params.currentUser}/${cdsID}`)
         .then(response => {
            this.currentCDS = response.data.cds;
            this.blastResults = response.data.blast;
            // let xmin = this.currentCDS.start-50;
            // let xmax = this.currentCDS.stop+50;
            this.startOptions = response.data.start_options;
            // this.series1 = [{ data: response.data.frame_1 }]
            // this.series2 = [{ data: response.data.frame_2 }]
            // this.series3 = [{ data: response.data.frame_3 }]
            // this.series4 = [{ data: response.data.frame_4 }]
            // this.series5 = [{ data: response.data.frame_5 }]
            // this.series6 = [{ data: response.data.frame_6 }]
            // this.options1 = {
            //    chart: { id: '1', group: 'direct', type: 'line', height: 170 },
            //    stroke: { width: 1.5 },
            //    colors: ['#008FFB'],
            //    noData: { text: 'Loading...' },
            //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
            //    xaxis: { min: xmin, max: xmax, tickAmount: 10 }
            // }
            // this.options2 = {
            //    chart: { id: '2', group: 'direct', type: 'line', height: 170 },
            //    stroke: { width: 1.5 },
            //    colors: ['#008FFB'],
            //    noData: { text: 'Loading...' },
            //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 },
            //             title: { text: 'Direct Sequences'}  },
            //    xaxis: { min: xmin, max: xmax, tickAmount: 10 }
            // }
            // this.options3 = {
            //    chart: { id: '3', group: 'direct', type: 'line', height: 193 },
            //    stroke: { width: 1.5 },
            //    colors: ['#008FFB'],
            //    noData: { text: 'Loading...' },
            //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
            //    xaxis: { min: xmin, max: xmax, tickAmount: 10,
            //             title: { text: 'Nucleotide Position'}  }
            // }
            // this.options4 = {
            //    chart: { id: '4', group: 'comp', type: 'line', height: 170 },
            //    stroke: { width: 1.5 },
            //    colors: ['#008FFB'],
            //    noData: { text: 'Loading...' },
            //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
            //    xaxis: { min: xmin, max: xmax, tickAmount: 10 }
            // }
            // this.options5 = {
            //    chart: { id: '5', group: 'comp', type: 'line', height: 170 },
            //    stroke: { width: 1.5 },
            //    colors: ['#008FFB'],
            //    noData: { text: 'Loading...' },
            //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 },
            //             title: { text: 'Complementary Sequences'} },
            //    xaxis: { min: xmin, max: xmax, tickAmount: 10 }
            // }
            // this.options6 = {
            //    chart: { id: '6', group: 'comp', type: 'line', height: 193 },
            //    stroke: { width: 1.5 },
            //    colors: ['#008FFB'],
            //    noData: { text: 'Loading...' },
            //    yaxis: { min: 0, max: 1, tickAmount: 5, labels: { minWidth: 40 } },
            //    xaxis: { min: xmin, max: xmax, tickAmount: 10,
            //             title: { text: 'Nucleotide Position'} }
            // }
         })
         .catch(error => {
            console.error(error);
         });
      },
      editCDS() {
         this.updatedCDS = this.currentCDS;
         this.updatedCDS.start = this.newStart;
         this.updatedCDS.function = this.newFunction;
         const payload = {
            id: this.updatedCDS.id,
            start: this.updatedCDS.start,
            stop: this.updatedCDS.stop,
            strand: this.updatedCDS.strand,
            function: this.updatedCDS.function,
            status: 'Pass',
         };
         this.updateCDS(payload, this.updatedCDS.id);
      },
      updateCDS(payload, cdsID) {
         axios.put(`http://localhost:5000/api/annotations/fail/${this.$route.params.currentUser}/${cdsID}`, payload)
         .then(() => {
            this.$router.push(`/annotations/${this.$route.params.currentUser}`);
         })
         .catch(error => {
            console.error(error);
         });
      },
      keepOriginal() {
         const payload = {
            id: this.currentCDS.id,
            start: this.currentCDS.start,
            stop: this.currentCDS.stop,
            strand: this.currentCDS.strand,
            function: 'None',
            status: 'Pass',
         };
         this.updateCDS(payload, this.currentCDS.id);
      },
      setFunction(funct) {
         this.newFunction = funct;
      },
      setStart(strt) {
         this.newStart = strt;
      },
   },
}
</script>

<style scoped>
/* ----- Title and Headers ----- */
.headers {
   margin: 40px auto;
}

.alert-primary {
   text-align: left;
   margin: 40px auto;
}

.subheader {
   text-align: left;
}

.info-bottom {
   margin: 50px auto;
}

/* ----- Column Styling ----- */
.columns {
   display: flex;
   flex-direction: row;
   justify-content: space-between;
}

.column {
   width: 100%;
}

/* ----- Coding Potential ----- */
#comp-graph {
   margin-bottom: 30px;
}

.table-responsive {
   max-height: 250px;
   overflow-y: auto;
   display: inline-block;
}

.table-responsive thead th {
   position: sticky;
   top: 0;
   background: #eee;
   border: darkgray;
}

caption {
  display: table-caption;
  caption-side: top;
}

tbody {
   width: 100%;
}

/* ----- Blast Results ----- */
.btn-blast {
   width: 100%;
}

/* --------------------- */
/* Mobile Styles */
@media only screen and (max-width: 1000px) {
   .columns {
      flex-direction: column;
   }

}
</style>