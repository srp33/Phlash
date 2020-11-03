<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :dnamaster="navDNAMaster"
      :blast="navBlast"
      :annotations="navAnnotations"
    />
    <div class="container">
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <h1>Annotations</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>
          Gene calls from DNA Master and GeneMark have been compared. A status
          has been assigned to each gene call from DNA Master.
        </p>
        <p>
          When the 'Action' column only contains 'Done' or 'Reinstate', you will
          be able to download your GenBank file.
        </p>
        <p><strong>Status Definitions</strong></p>
        <p>
          <strong style="color: green">Green:</strong> Pass, the open reading frame is
          the same as or longer than GeneMark's gene call.<br />
          <strong style="color: red">Red:</strong> Fail, the open reading frame is
          shorter than GeneMark's gene call.<br />
          <strong style="color: orange">Orange:</strong> Undetermined, the open reading frame
          has not been called at all in GeneMark or the stop position differs from GeneMark's gene call.<br />
          <strong>S: </strong>Short, the open reading frame is less than 
          200 base pairs.<br />
          <strong>LLG: </strong>Long leading gap, there is a large gap between this 
          and the previous gene<br />
          <strong>LLO: </strong>Long leading overlap, the gene overlaps the previous 
          gene by more than 10 base pairs.<br />
          <strong>LTG: </strong>Long tailing gap, the gene overlaps the next 
          gene by more than 10 base pairs.<br />
          <strong>LTO: </strong>Long tailing overlap, the gene overlaps the next  
          gene by more than 10 base pairs.<br />
          <strong>STG: </strong>Short tailing gap, the gene is on a different strand 
          than the previous gene and has a gap less than 50 base pairs in length.<br />
          <strong>SLG: </strong>Short leading gap, the gene is on a different strand 
          than the next gene and has a gap less than 50 base pairs in length.<br />
          <strong>GOOD: </strong>The gene has a reasonable length and overlap.     
        </p>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Blast', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-light btn-nav">
              <svg
                class="bi bi-arrow-left"
                width="1em"
                height="1em"
                viewBox="0 0 16 16"
                fill="currentColor"
                xmlns="http://www.w3.org/2000/svg"
              >
                <path
                  fill-rule="evenodd"
                  d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z"
                  clip-rule="evenodd"
                />
                <path
                  fill-rule="evenodd"
                  d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z"
                  clip-rule="evenodd"
                />
              </svg>
              <strong>Back</strong>
            </button>
          </router-link>
        </div>
      </div>
      <div class="alert alert-dark" v-if="completedGenes != dnamaster.length">
        You have
        <strong
          >{{ dnamaster.length - completedGenes }}/{{
            dnamaster.length
          }}</strong
        >
        genes remaining.
      </div>
      <div
        class="alert alert-success"
        v-if="completedGenes == dnamaster.length"
      >
        Congratulations! You can now
        <a href="#" @click="downloadGenBankFile" class="alert-link">
          download your GenBank file</a
        >.
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
                <td v-else-if="curr.function.length < 16">
                  {{ curr.function }}
                </td>
                <td v-else>{{ curr.function.substring(0, 13) }}...</td>
                <td v-if="curr.function == 'DELETED'"></td>
                <td v-else-if="curr.status == 'Pass'" style="color: green">
                  {{ getStatus(index) }}
                </td>
                <td v-else-if="curr.status == 'Fail'" style="color: red">
                  {{ getStatus(index) }}
                </td>
                <td v-else style="color: orange">
                  {{ getStatus(index) }}
                </td>
                <td v-if="curr.function == 'DELETED'">
                  <button
                    class="btn btn-outline-dark btn-sm"
                    style="width: 100px"
                    @click="reinstate(index)"
                  >
                    <strong>Reinstate</strong>
                  </button>
                </td>
                <td v-else>
                  <router-link
                    :to="{
                      name: 'CDS',
                      params: {
                        phageID: $route.params.phageID,
                        cdsID: curr.id,
                      },
                    }"
                  >
                    <button
                      class="btn btn-outline-dark btn-sm"
                      style="width: 100px"
                      v-if="curr.function !== 'None selected'"
                    >
                      <strong>Done</strong>
                    </button>
                  </router-link>
                  <router-link
                    :to="{
                      name: 'CDS',
                      params: {
                        phageID: $route.params.phageID,
                        cdsID: curr.id,
                      },
                    }"
                  >
                    <button
                      class="btn btn-dark btn-sm"
                      style="width: 100px"
                      v-if="curr.function === 'None selected'"
                    >
                      <strong>Go</strong>
                    </button>
                  </router-link>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
      <div class="nav-btns-wrapper">
        <router-link
          :to="{ name: 'Blast', params: { phageID: $route.params.phageID } }"
        >
          <button class="btn btn-light btn-nav">
            <svg
              class="bi bi-arrow-left"
              width="1em"
              height="1em"
              viewBox="0 0 16 16"
              fill="currentColor"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                fill-rule="evenodd"
                d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z"
                clip-rule="evenodd"
              />
              <path
                fill-rule="evenodd"
                d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z"
                clip-rule="evenodd"
              />
            </svg>
            <strong>Back</strong>
          </button>
        </router-link>
      </div>
      <b-modal
        v-model="showModal"
        id="my-modal"
        title="Congratulations!"
        hide-footer
      >
        <p>You have finished your phage genome annotations!</p>
        <b-button
          class="mt-3"
          block
          style="margin-top: 0px"
          @click="downloadGenBankFile"
        >
          <strong>Download GenBank file</strong>
          <div
            v-if="downloadLoading"
            class="spinner-border spinner-border-sm"
            role="status"
          >
            <span class="sr-only"></span>
          </div>
        </b-button>
      </b-modal>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "Annotations",
  components: {
    Loading,
    Navbar,
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
        strand: "",
      },
      fileDownloaded: false,
      pageLoading: true,
      downloadLoading: false,
      completedGenes: 0,
    };
  },

  created() {
    this.getData();
  },

  computed: {

    navUpload: function () {
      return true;
    },

    navDNAMaster: function () {
      return true;
    },

    navBlast: function () {
      return true;
    },

    navAnnotations: function () {
      return true;
    },

  },

  methods: {

    /**
     * Populates dnamaster with data.
     */
    getData() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}`
        )
        .then((response) => {
          this.dnamaster = response.data.dnamaster;
          this.pageLoading = false;
          for (var i = 0; i < this.dnamaster.length; i++) {
            if (this.dnamaster[i].function != "None selected")
              ++this.completedGenes;
          }
          if (this.completedGenes == this.dnamaster.length) {
            this.showModal = true;
          }
        })
        .catch((error) => {
          console.error(error);
        });
    },
    
    /**
     * Downloads the GenBank file upon annotation completion.
     */
    downloadGenBankFile() {
      this.downloadLoading = true;
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}`
        )
        .then((response) => {
          let data = response.data;
          const blob = new Blob([data]);
          let link = document.createElement("a");
          link.href = window.URL.createObjectURL(blob);
          link.download = `${this.$route.params.phageID}_modified.gb`;
          link.click();
          this.downloadLoading = false;
          this.fileDownloaded = true;
        });
    },

    /**
     * Changes the function of a deleted gene so that it is now visible.
     */
    reinstate(index) {
      this.completedGenes -= 1;
      this.dnamaster[index].function = "None selected";
    },

    /**
     * Creates the status for a given cds.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getStatus(index) {
      var status = "";
      if (this.dnamaster[index].stop - this.dnamaster[index].start < 200) {
        status += "S | ";
      }
      if (index == 0) {
        status += this.getTailingOverlap(index);
      }
      else if (index == this.dnamaster.length - 1) {
        status += this.getLeadingOverlap(index);
      }
      else {
        status += this.getLeadingOverlap(index);
        status += this.getTailingOverlap(index);
      }
      if (status == "") { return "GOOD"; }
      return status.substring(0, status.length - 3);
    },

    /**
     * Determines if there is a large gap or overlap between the current gene and the next.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getTailingOverlap(index) {
      var status = "";
      var overlap = this.dnamaster[index].stop - this.dnamaster[index + 1].start
      if (this.dnamaster[index].strand == this.dnamaster[index + 1].strand) {
        if (overlap < -10) {
          status += "LTG | ";
        }
        else if (overlap > 10) {
          status += "LTO | ";
        }
      }
      else if (overlap > -50) {
        status += "STG | "
      }
      return status;
    },

    /**
     * Determines if there is a large gap or overlap between the current gene and the previous.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getLeadingOverlap(index) {
      var status = "";
      var leadingOverlap = this.dnamaster[index - 1].stop - this.dnamaster[index].start
      if (this.dnamaster[index].strand == this.dnamaster[index - 1].strand) {
        if (leadingOverlap < -10) {
          status += "LLG | ";
        }
        else if (leadingOverlap > 10) {
          status += "LLO | ";
        }
      }
      else if (leadingOverlap > -50) {
      status += "SLG | "
      }
      return status;
    }  
  },
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
