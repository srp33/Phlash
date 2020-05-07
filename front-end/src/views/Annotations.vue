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
          has been assigned to each gene call from DNA Master. To understand what
          each status means, please refer to the definitions below.
        </p>
        <p>
          When the 'Action' column only contains 'Done', you will be able to
          download your GenBank file.
        </p>
        <button class="btn btn-light" style="margin-top: 0px;" @click="downloadGenBankFile">
          <strong>Download GenBank file</strong>
            <div v-if="downloadLoading" class="spinner-border spinner-border-sm" role="status">
              <span class="sr-only"></span>
            </div>
        </button>
        <p><strong>Status Definitions</strong></p>
        <p>
          <strong style="color:green">Pass:</strong> DNA Master's gene call is the same as or longer than GeneMark's gene call.<br />
          <strong style="color:red">Fail:</strong> DNA Master's gene call is shorter than GeneMark's gene call.<br />
          <strong style="color:orange">Need more information:</strong> DNA Master's gene call has not been called at all in GeneMark.<br />
        </p>
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
                <td>{{ curr.id }}</td>
                <td>{{ curr.start }}</td>
                <td>{{ curr.stop }}</td>
                <td>{{ curr.strand }}</td>
                <td v-if="curr.function.length<30">{{ curr.function }}</td>
                <td v-else>{{ curr.function.substring(0,30) }}...</td>
                <td v-if="curr.status == 'Pass'" style="color: green">{{ curr.status }}</td>
                <td v-if="curr.status == 'Fail'" style="color: red">{{ curr.status }}</td>
                <td v-if="curr.status == 'Need more information'" style="color: orange">
                  {{ curr.status }}
                </td>
                <td>
                  <button class="btn btn-dark btn-sm" disabled style="width:100px"
                    v-if="curr.function !== 'None'">
                    <strong>Done</strong>
                  </button>
                  <router-link :to="{ name: 'CDS', params: {phageID: $route.params.phageID, cdsID: curr.id} }">
                    <button class="btn btn-dark btn-sm" style="width:100px"
                      v-if="curr.function === 'None'">
                      <strong>Go</strong>
                    </button>
                  </router-link>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
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
      axios.get(`http://localhost:5000/api/annotations/${this.$route.params.phageID}`)
        .then(response => {
          this.dnamaster = response.data.dnamaster;
          this.pageLoading = false;
        })
        .catch(error => {
          console.error(error);
        });
    },
    downloadGenBankFile() {
      this.downloadLoading = true;
      axios.post(`http://localhost:5000/api/annotations/${this.$route.params.phageID}`)
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

.btn-light {
  margin: 15px auto;
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