<template>
  <div class="wrapper">
    <Navbar :upload="navUpload" :dnamaster="navDNAMaster" :blast="navBlast" :annotations="navAnnotations" />
    <div class="container">
      <loading :active.sync="pageLoading" :is-full-page="true" :height="100" :width="100"></loading>
      <h1>Annotate</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>
          Gene calls from DNA Master and GeneMark have been compared. A status
          has been assigned to each gene call from DNA Master.
        </p>
        <p>
          When the 'Action' column only contains 'Done', you will be able to
          download your GenBank file.
        </p>
        <p><strong>Status Definitions</strong></p>
        <p>
          <strong style="color:green">Pass:</strong> DNA Master's gene call is the same as or longer than GeneMark's gene call.<br />
          <strong style="color:red">Fail:</strong> DNA Master's gene call is shorter than GeneMark's gene call.<br />
          <strong style="color:orange">Undetermined:</strong> DNA Master's gene call has not been called at all in GeneMark.<br />
        </p>
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Blast', params: {phageID: $route.params.phageID} }">
            <button class="btn btn-light btn-nav">
              <svg class="bi bi-arrow-left" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
                <path fill-rule="evenodd" d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z" clip-rule="evenodd"/>
                <path fill-rule="evenodd" d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z" clip-rule="evenodd"/>
              </svg>
              <strong>Back</strong>
            </button>
          </router-link>
        </div>
      </div>
      <div class="alert alert-warning" v-if="completed != dnamaster.length">
        You have <strong>{{ dnamaster.length - completed }}/{{ dnamaster.length }}</strong> genes remaining.
      </div>
      <div class="alert alert-success" v-if="completed == dnamaster.length">
        <strong> Congratulations! You can now download your GenBank file below. </strong>
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
                <th width="200px" scope="col">Function</th>
                <th scope="col">Status</th>
                <th scope="col">Action</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(curr, index) in dnamaster" :key="index">
                  <td v-if="curr.function == 'DELETED'">{{ curr.id }}</td>
                  <td v-else>{{ curr.id }}</td>
                  <td v-if="curr.function == 'DELETED'"></td>
                  <td v-else>{{ curr.start }}</td>
                  <td v-if="curr.function == 'DELETED'"></td>
                  <td v-else>{{ curr.stop }}</td>
                  <td v-if="curr.function == 'DELETED'"></td>
                  <td v-else>{{ curr.strand }}</td>
                  <td v-if="curr.function == 'DELETED'"></td>
                  <td v-else-if="curr.function.length<30">{{ curr.function }}</td>
                  <td v-else>{{ curr.function.substring(0,30) }}...</td>
                  <td v-if="curr.function == 'DELETED'"></td>
                  <td v-else-if="curr.status == 'Pass'" style="color: green">{{ curr.status }}</td>
                  <td v-else-if="curr.status == 'Fail'" style="color: red">{{ curr.status }}</td>
                  <td v-else-if="curr.status == 'Undetermined'" style="color: orange">{{ curr.status }}</td>
                  <td v-if="curr.function == 'DELETED'">
                    <router-link :to="{ name: 'CDS', params: {phageID: $route.params.phageID, cdsID: curr.id} }">
                      <button class="btn btn-dark btn-sm" style="width:100px">
                        <strong>Reinstate</strong>
                      </button>
                    </router-link>
                  </td>
                  <td v-else>
                    <router-link :to="{ name: 'CDS', params: {phageID: $route.params.phageID, cdsID: curr.id} }">
                      <button class="btn btn-dark btn-sm" disabled style="width:100px"
                        v-if="curr.function !== 'None selected'">
                        <strong>Done</strong>
                      </button>
                    <router-link :to="{ name: 'CDS', params: {phageID: $route.params.phageID, cdsID: curr.id} }">
                      <button class="btn btn-dark btn-sm" style="width:100px"
                        v-if="curr.function === 'None selected'">
                        <strong>Go</strong>
                      </button>
                    </router-link>
                  </td> 
              </tr>
            </tbody>
          </table>
        </div>
      </div>
      <button class="btn btn-light btn-gb" v-if="completed == dnamaster.length" style="margin-top: 0px;" @click="downloadGenBankFile">
        <strong>Download GenBank file</strong>
          <div v-if="downloadLoading" class="spinner-border spinner-border-sm" role="status">
            <span class="sr-only"></span>
          </div>
      </button>
      <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Blast', params: {phageID: $route.params.phageID} }">
            <button class="btn btn-light btn-nav">
              <svg class="bi bi-arrow-left" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
                <path fill-rule="evenodd" d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z" clip-rule="evenodd"/>
                <path fill-rule="evenodd" d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z" clip-rule="evenodd"/>
              </svg>
              <strong>Back</strong>
            </button>
          </router-link>
        </div>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue"
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "Annotations",
  components: {
    Loading,
    Navbar
  },
  data() {
    return {
      dnamaster: [],
      startOptions: [],
      showModal: false,
      currCDS: {
        id: "",
        start: "",
        stop: "",
        strand: ""
      },
      fileDownloaded: false,
      pageLoading: true,
      downloadLoading: false,
      completed: 0,
    };
  },
  created() {
    this.getData();
  },
  computed: {
    navUpload: function() {
      return true;
    },
    navDNAMaster: function() {
      return true;
    },
    navBlast: function() {
      return true;
    },
    navAnnotations: function() {
      return true;
    },
  },
  methods: {
    getData() {
      axios.get(process.env.VUE_APP_BASE_URL + `/annotations/${this.$route.params.phageID}`)
        .then(response => {
          this.dnamaster = response.data.dnamaster;
          this.pageLoading = false;
          for (var i = 0; i < this.dnamaster.length; i++) {
            if (this.dnamaster[i].function != "None selected") ++this.completed;
          }
        })
        .catch(error => {
          console.error(error);
        });
    },
    downloadGenBankFile() {
      this.downloadLoading = true;
      for (var i = 0; i < this.dnamaster.length; i++) {
            if (this.dnamaster[i].function == "DELETED") {
              axios.delete(process.env.VUE_APP_BASE_URL + `/annotations/cds/${this.$route.params.phageID}/${this.dnamaster[i].id}`)
                .then(() => {
                  this.$router.push(`/annotations/${this.$route.params.phageID}`);
                })
                .catch(error => {
                  console.error(error);
                });
            }
      }
      axios.post(process.env.VUE_APP_BASE_URL + `/annotations/${this.$route.params.phageID}`)
        .then(response => {
          let data = response.data;
          const blob = new Blob([data]);
          let link = document.createElement("a");
          link.href = window.URL.createObjectURL(blob);
          link.download = `${this.$route.params.phageID}_modified.gb`;
          link.click();
          this.downloadLoading = false;
          this.fileDownloaded = true;
        });
    }
  }
};
</script>

<style scoped>
.wrapper {
  margin: 0;
}

h1 {
  margin: 40px auto;
}

.alert-primary {
  text-align: left;
  margin: 40px auto;
}

.btn-gb {
  margin: 15px auto;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 10px;
}

.bi-arrow-left {
  margin-right: 5px;
  margin-left: 0px;
}

.status-btn {
  width: 100%;
}

td,
th {
  word-wrap: break-word;
  width: 150px;
}
</style>
