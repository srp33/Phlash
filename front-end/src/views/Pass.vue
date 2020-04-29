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
         <p>Choose a function for this gene call from the BLAST results below.</p>
         <p><strong>Your function selection:</strong> {{ newFunction }}</p>
         <button type="button" class="btn btn-light" @click="editCDS" v-if="newFunction!=='None'">
            <strong>Update {{ $route.params.id }}</strong>
         </button>
         <button type="button" class="btn btn-light" @click="deleteCDS($route.params.id)" v-if="!blastResultsExist">
            <strong>Delete {{ $route.params.id }}</strong>
         </button>
      </div>
      <div class="blast-results">
         <h3 style="text-align: center; margin: 20px;">BLAST Results</h3>
         <BlastResults :blastResults="blastResults" :blastResultsExist="blastResultsExist" @newFunction="setFunction" />
      </div>
      <div class="info-bottom">
         <p><strong>Your function selection:</strong> {{ newFunction }}</p>
         <button type="button" class="btn btn-light" @click="editCDS" v-if="newFunction!=='None'">
            <strong>Update {{ $route.params.id }}</strong>
         </button>
         <button type="button" class="btn btn-light" @click="deleteCDS($route.params.id)" v-if="!blastResultsExist">
            <strong>Delete {{ $route.params.id }}</strong>
         </button>
      </div>
   </div>
</template>

<script>
import axios from 'axios';
import BlastResults from '../components/BlastResults.vue'

export default {
   name: 'Pass',
   components: {
      BlastResults
   },
   data () {
      return {
         blastResults: [],
         blastResultsExist: true,
         showFunction: false,
         newFunction: 'None',
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
      }
   },
   created() {
      this.getData(this.$route.params.id);
   },
   methods: {
      getData(cdsID) {
         axios.get(`http://localhost:5000/api/annotations/pass/${this.$route.params.currentUser}/${cdsID}`)
         .then(response => {
            this.currentCDS = response.data.cds;
            this.blastResults = response.data.blast;
            if (this.blastResults.length === 0) this.blastResultsExist = false;
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
         axios.put(`http://localhost:5000/api/annotations/pass/${this.$route.params.currentUser}/${cdsID}`, payload)
         .then(() => {
            this.$router.push(`/annotations/${this.$route.params.currentUser}`);
         })
         .catch(error => {
            console.error(error);
         });
      },
      deleteCDS(cdsID) {
         axios.delete(`http://localhost:5000/api/annotations/pass/${this.$route.params.currentUser}/${cdsID}`)
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
            function: this.currentCDS.function,
            status: this.currentCDS.status,
         };
         this.updateCDS(payload, this.currentCDS.id);
      },
      checkBlastResults() {
         if (this.blastResults) {
            if (this.blastResults.length > 0) return true
            else return false
         } else {
            return false
         }
      }
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

/* ----- ----- */
.info-bottom {
   margin: 50px auto;
}
</style>