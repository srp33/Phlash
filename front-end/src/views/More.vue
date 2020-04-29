<template>
<div class="container">
   <loading loader="dots" :active.sync="isLoading" :is-full-page="true" :height=200 :width=200 color="#2279b6"></loading>
   <div class="headers">
      <h1>ID: {{ $route.params.id }}</h1>
      <h4>Start: {{ currentCDS.start }}</h4>
      <h4>Stop: {{ currentCDS.stop }}</h4>
      <h4>Status: {{ currentCDS.status }}</h4>
   </div>
   <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>
         Based on the information given below, you may choose to keep this gene
         call and choose a new function for it, or delete the gene call altogether.
      </p>
      <p><strong>Your selected function:</strong> {{ newFunction }}</p>
      <button type="button" class="btn btn-light btn-action" @click="editCDS" v-if="newFunction !== 'None'">
         <strong>Update {{ $route.params.id }}</strong>
      </button>
      <button type="button" class="btn btn-light btn-action" @click="deleteCDS($route.params.id)">
         <strong>Delete {{ $route.params.id }}</strong>
      </button>
   </div>
   <div class="columns-wrapper">
      <h3 style="text-align: center; margin: 40px;">GeneMark's Coding Potential</h3>
      <div class="subheader">

      </div>
   </div>
   <div class="coding-potential-table">
      <h5 style="text-align: center; margin: 20px;">Average Coding Potential Per Frame</h5>
      <div class="table-responsive">
         <table id="cp-table" class="table table-hover">
            <caption>Frames are counted 1-6 (direct 1-3 and complementary 4-6).</caption>
            <thead>
               <tr>
                  <th scope="col">Start Position</th>
                  <th scope="col">Frame 1</th>
                  <th scope="col">Frame 2</th>
                  <th scope="col">Frame 3</th>
                  <th scope="col">Frame 4</th>
                  <th scope="col">Frame 5</th>
                  <th scope="col">Frame 6</th>
               </tr>
            </thead>
            <tbody>
               <tr>
                  <th>{{ currentCDS.start }}</th>
                  <td>{{ probabilities.frame_1 }}</td>
                  <td>{{ probabilities.frame_2 }}</td>
                  <td>{{ probabilities.frame_3 }}</td>
                  <td>{{ probabilities.frame_4 }}</td>
                  <td>{{ probabilities.frame_5 }}</td>
                  <td>{{ probabilities.frame_6 }}</td>
               </tr>
            </tbody>
         </table>
      </div>
   </div>
   <div class="coding-potential-graphs">
      <h5 style="text-align: center; margin: 20px;">Coding Potential Per Frame</h5>
      <div v-if="dataExists">
         <Graphs :data1="data1" :data2="data2" :data3="data3" :data4="data4" :data5="data5" :data6="data6" :start="currentCDS.start" :stop="currentCDS.stop"/>
      </div>
      <div v-else>
         <p>Loading graphs...</p>
      </div>
   </div>
   <hr/>
   <div class="blast-results">
      <h3 style="text-align: center; margin: 20px;">BLAST Results</h3>
      <BlastResults :blastResults="blastResults" :blastResultsExist="blastResultsExist" @newFunction="setFunction" />
   </div>
   <div class="info-bottom">
      <p><strong>Your function selection:</strong> {{ newFunction }}</p>
      <button type="button" class="btn btn-light btn-action" @click="editCDS" v-if="newFunction !== 'None'">
         <strong>Update {{ $route.params.id }}</strong>
      </button>
      <button type="button" class="btn btn-light btn-action" @click="deleteCDS($route.params.id)">
         <strong>Delete {{ $route.params.id }}</strong>
      </button>
   </div>
</div>
</template>

<script>
import BlastResults from '../components/BlastResults.vue'
import Graphs from '../components/Graphs.vue'
import axios from 'axios'
import Loading from 'vue-loading-overlay'
import 'vue-loading-overlay/dist/vue-loading.css';

export default {
   name: 'More',
   components: {
      Loading,
      Graphs,
      BlastResults
   },
   data() {
      return {
         probabilities: [],
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
         blastResultsExist: true,
         newFunction: 'None',
         dataExists: false,
         isLoading: true,
         data1: [{ x:[], y:[] }],
         data2: [{ x:[], y:[] }],
         data3: [{ x:[], y:[] }],
         data4: [{ x:[], y:[] }],
         data5: [{ x:[], y:[] }],
         data6: [{ x:[], y:[] }],
      }
   },
   created() {
      this.getData(this.$route.params.id);
   },
    methods: {
      getData(cdsID) {
         axios.get(`http://localhost:5000/api/annotations/more/${this.$route.params.currentUser}/${cdsID}`)
         .then(response => {
            this.currentCDS = response.data.cds;
            this.blastResults = response.data.blast;
            if (this.blastResults.length === 0) this.blastResultsExist = false;
            this.probabilities = response.data.probs;
            this.data1 = [{ 
               x: response.data.x_data,
               y: response.data.y_data_1
            }]
            this.data2 = [{ 
               x: response.data.x_data,
               y: response.data.y_data_2
            }]
            this.data3 = [{ 
               x: response.data.x_data,
               y: response.data.y_data_3
            }]
            this.data4 = [{ 
               x: response.data.x_data,
               y: response.data.y_data_4
            }]
            this.data5 = [{ 
               x: response.data.x_data,
               y: response.data.y_data_5
            }]
            this.data6 = [{ 
               x: response.data.x_data,
               y: response.data.y_data_6
            }]
            this.dataExists = true;
            this.isLoading = false;
         })
         .catch(error => {
            console.error(error);
         });
      },
      setFunction(funct) {
         this.newFunction = funct;
      },
      editCDS() {
         this.updatedCDS = this.currentCDS;
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
         axios.put(`http://localhost:5000/api/annotations/more/${this.$route.params.currentUser}/${cdsID}`, payload)
         .then(() => {
            this.$router.push(`/annotations/${this.$route.params.currentUser}`);
         })
         .catch(error => {
            console.error(error);
         });
      },
      deleteCDS(cdsID) {
         axios.delete(`http://localhost:5000/api/annotations/more/${this.$route.params.currentUser}/${cdsID}`)
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

.btn-action {
   margin: auto 10px;
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
  /* caption-side: top; */
}

tbody {
   width: 100%;
}

.coding-potential-graphs {
   height: 670px;
}

/* ----- Blast Results ----- */
.btn-blast {
   width: 100%;
}

/* --------------------- */
/* Mobile Styles */
@media only screen and (max-width: 1200px) {
   .coding-potential-graphs {
      height: 1300px;
   }
}
</style>