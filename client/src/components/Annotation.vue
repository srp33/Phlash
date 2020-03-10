<template>
<div class="container">
   <h1><strong>Annotate:</strong> Blast</h1>
   <h3>DNA Master vs. GeneMark</h3>
   <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>
         <strong style="color:green">Pass:</strong> DNA Master's gene call is the same as or longer than GeneMark's gene call.<br />
         <strong style="color:red">Fail:</strong> DNA Master's gene call is shorter than GeneMark's gene call.<br /> 
         <strong style="color:orange">Need more information:</strong> DNA Master's gene call has not been called at all in GeneMark.<br/>
      </p>
      <p>For gene calls that have the `Fail` status, click on each button and choose a start position to use for BLAST.</p>
      <button class="btn btn-light" @click="downloadFile">
         <strong>Download FASTA file for BLAST</strong>
      </button>
      <button class="btn btn-light">
         <strong>Go to NCBI's BLASTp</strong>
      </button>
   </div>
   <div id="annotations" align="center">
      <div class="table-responsive">
         <table class="table table-hover" align="center">
               <thead>
                  <tr>
                     <th scope="col">ID</th>
                     <th scope="col">Start</th>
                     <th scope="col">Stop</th>
                     <th scope="col">Strand</th>
                     <th scope="col">Function</th>
                     <th scope="col">Status</th>
                     <!-- <th scope="col">Action</th> -->
                  </tr>
               </thead>
               <tbody>
                  <tr v-for="(curr, index) in dnamaster" :key="index">
                     <td>{{ curr.id }}</td>
                     <td>{{ curr.start }}</td>
                     <td>{{ curr.stop }}</td>
                     <td>{{ curr.strand }}</td>
                     <td v-if="curr.function.length<30">{{ curr.function }}</td>
                     <td v-else>{{ curr.function.substring(0,30) }}...</td>
                     <td v-if="curr.status == 'Pass'" style="color: green">{{ curr.status }}</td>
                     <td v-if="curr.status == 'Fail'" style="color: red">{{ curr.status }}</td>
                     <td v-if="curr.status == 'Need more information'" style="color: orange">{{ curr.status }}</td>
                     <!-- <td>{{ curr.status }}</td> -->
                     <!-- <td>
                        <span style="color:green" 
                              v-if="curr.function !== 'None' && curr.status == 'Pass'">
                           <strong>{{ curr.status }}</strong>
                        </span>
                        <router-link :to="{ name: 'Blast', params: {id: curr.id} }">
                           <button class="btn btn-success btn-sm" 
                                    v-if="curr.function == 'None' && curr.status == 'Pass'">
                              <strong>{{ curr.status }}</strong>
                           </button>
                        </router-link>
                        <router-link :to="{ name: 'Failed', params: {id: curr.id} }">
                           <button class="btn btn-danger btn-sm" style="width:141px"
                                    v-if="curr.status == 'Fail'">
                              <strong>{{ curr.status }}</strong>
                           </button>
                        </router-link>
                        <router-link :to="{ name: 'More', params: {id: curr.id} }">
                           <button class="btn btn-warning btn-sm" style="width:141px"
                                    v-if="curr.status =='Need more information'">
                              <strong>{{ curr.status }}</strong>
                           </button>
                        </router-link>
                     </td> -->
                  </tr>
               </tbody>
         </table>
      </div>
   </div>
   <b-modal ref="chooseStartModal" id="modal" v-model="showModal" title="currCDS.id: currCDS.start - currCDS.stop">
      <b-form @submit="onSubmit" @reset="onReset" align="left">
         <p>Choose additional start positions to BLAST:</p>
         <div v-for="start in startOptions" :key="start">
            <div class="form-check">
               <input class="form-check-input" type="checkbox" value="">
               <label class="form-check-label">
                  {{ start.start_position }} - {{ curr.CDS.stop }}
               </label>
            </div>
         </div>
      </b-form>
   </b-modal>
</div>
</template>

<script>
import axios from 'axios';

export default {
   data() {
      return {
         dnamaster: [],
         showFail: false,
         startOptions: [],
         showModal: false,
         currCDS: {
            id: '',
            start: '',
            stop: '',
            strand: '',
         },
      };
   },
   methods: {
      getData() {
         axios.get('http://localhost:5000/annotate')
         .then(response => {
            this.dnamaster = response.data.dnamaster;
         })
         .catch(error => {
            console.error(error);
         });
      },
      downloadFile() {
         axios.post('http://localhost:5000/annotate')
         .then(response => {
            let data = response.data;
            const blob = new Blob([data], { type: 'application/fasta' })
            let link = document.createElement('a')
            link.href = window.URL.createObjectURL(blob)
            link.download = 'sequences_to_blast.fasta'
            link.click()
         });
      },
      getStartData(cds) {
         this.currCDS.id = cds.id;
         this.currCDS.start = cds.start;
         this.currCDS.stop = cds.stop;
         this.currCDS.strand = cds.strand;

         axios.get(`http://localhost:5000/annotate/failed/${cds.id}`)
         .then(response => {
            this.startOptions = response.data.start_options;
            this.showModal = true;
         })
         .catch(error => {
            console.error(error);
         });
      },
      setCurrCDS(cds) {
         this.currCDS.id = cds.id;
         this.currCDS.start = cds.start;
         this.currCDS.stop = cds.stop;
         this.currCDS.strand = cds.strand;
         this.getStartData(this.currCDS.id);
      },
      initForm() {
         this.currCDS.id = '';
         this.currCDS.start = '';
         this.currCDS.stop = '';
         this.currCDS.strand = '';
      },
      onSubmit(evt) {
         evt.preventDefault();
         this.$refs.chooseStartModal.hide();
         // FIXME: ADD CHOSEN STARTS TO BLAST FASTA FILE
         this.initForm()
      },
      onReset(evt) {
         evt.preventDefault();
         this.$refs.chooseStartModal.hide();
         this.initForm();
         this.getData();
      },
   },
   created() {
      this.getData();
   },
};
</script>

<style scoped>
h2 {
   margin-bottom: 40px;
}

.alert-primary {
   text-align: left;
   margin: 40px auto;
}

.btn-light {
   margin: 5px;
}

.status-btn {
   width: 100%;
}

td, th {
   word-wrap: break-word;
   width: 150px;
}
</style>