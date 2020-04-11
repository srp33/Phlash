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
      <h3 style="text-align: center; margin: 20px;">Coding Potential</h3>
      <div class="subheader">

      </div>
   </div>
   <div class="coding-potential-table">
      <div class="subheader">
         <p>
            Average coding potential per
            frame between the specified start position and {{
            currentCDS.stop }}. Frames are counted 1-6
            (direct 1-3 and complementary 4-6).
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
   <hr/>
   <div class="blast-results">
      <h3 style="text-align: center; margin: 20px;">BLAST Results</h3>
      <BlastResults :blastResults="blastResults" @newFunction="setFunction" />
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
import axios from 'axios'
// import VueApexCharts from 'vue-apexcharts'

export default {
   name: 'More',
   components: {
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
         newFunction: 'None',
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
            this.probabilities = response.data.probs;
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