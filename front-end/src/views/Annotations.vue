<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
      :logout="true"
    />
    <div class="container">
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <h1>Annotations</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <p>
          Gene calls that have been made from <i>GeneMarkS</i>, <i>Glimmer3</i>,
          <i>ARAGORN</i>, and <i>PHANOTATE</i> are shown below. To have the most
          accurate annotation of this genome manual annotation is necessary. The
          status for each gene is shown below which are mearly guides and
          suggestions. The status does not actually mean that a gene call is
          correct or incorrect. The parameters used for creating the status for
          each gene can be updated in the 'settings' tab above.
        </p>
        <hr />
        <p><strong>Status Definitions</strong></p>
        <p>
          <strong style="color: #1b9e77">Green:</strong> Pass, the open reading
          frame covers the coding potential.<br />
          <strong style="color: #d95f02">Orange:</strong> Fail, the open reading
          frame does not cover the coding potential.<br />
          <strong style="color: #7570b3">Purple:</strong> tRNA, this gene
          represents a tRNA.<br />
          <strong>S: </strong>Short, the open reading frame is less than
          {{ short }} base pairs.<br />
          <strong>LLG: </strong>Long leading gap, there is more than a
          {{ gap }} base pair gap between this and the previous gene.<br />
          <strong>LLO: </strong>Long leading overlap, the gene overlaps the
          previous gene by more than {{ overlap }} base pairs.<br />
          <strong>LTG: </strong>Long tailing gap, there is more than a
          {{ gap }} base pair gap between this and the next gene.<br />
          <strong>LTO: </strong>Long tailing overlap, the gene overlaps the next
          gene by more than {{ overlap }} base pairs.<br />
          <strong>STG: </strong>Short tailing gap, the gene is on a different
          strand than the previous gene and has a gap less than
          {{ oppositeGap }} base pairs in length.<br />
          <strong>SLG: </strong>Short leading gap, the gene is on a different
          strand than the next gene and has a gap less than
          {{ oppositeGap }} base pairs in length.<br />
          <strong>R: </strong>The gene has a reasonable length and overlap.
        </p>
        <hr />
        <p><strong>Actions</strong></p>
        <p v-if="!viewOnly">
          <strong>Annotate: </strong>When clicked, a function and alternate open
          reading frame can be added.<br />
          <strong>Delete: </strong>When clicked, this gene will be temporarily
          removed.<br />
          <strong>Reinstate: </strong>When clicked, the deleted gene will be
          added.<br />
          <strong>White background: </strong>This CDS has already been updated
          but can still be edited.<br />
          In rare occasions when a needed coding sequence is not shown, the 'Add
          CDS' button may be clicked to add a new custom CDS.<br />
          <button v-if="!viewOnly" class="btn btn-dark btn-action" @click="showAddCDS = true">
            <strong>&#43; Add CDS</strong>
          </button>
        </p>
        <p v-else>
          <strong>View: </strong>When clicked, the annotations for this CDS may be viewed.<br />
          <strong>White background: </strong>This CDS has already been updated by the owner.
        </p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Blast', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'GeneMap',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav" id="next-top">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
      <div v-if="blastLoading" class="alert alert-primary alert-dismissible text-size">
        <a href="#" class="close" data-dismiss="alert" aria-label="close"
          >&times;</a
        >
        {{waitMessage}}<br />
        You may continue to work on this phage's genome, but keep in mind that all the Blast results will not be shown. 
        If 'Annotate' is clicked for a CDS that does not currently have any
        data, you will be brought back to this page.<br />
        A notification will appear when finished.
      </div>
      <div
        class="alert alert-secondary"
        style="text-align: center"
        v-if="completedGenes !== phageAnnotations.length"
      >
        You have
        <strong
          >{{ phageAnnotations.length - completedGenes }}/{{
            phageAnnotations.length
          }}</strong
        >
        genes remaining.
      </div>
      <div
        class="alert alert-secondary"
        v-if="completedGenes === phageAnnotations.length"
      >
        Congratulations! You have annotated every CDS. Click 'Next' to see a map
        of the genome.
      </div>
      <div id="annotations" align="center">
        <div
          class="table table-responsive table-secondary"
          style="overflow-y: auto; max-height: 50em"
        >
          <table class="table table-hover" align="center">
            <thead>
              <tr>
                <th scope="col">ID</th>
                <th scope="col">Left</th>
                <th scope="col">Right</th>
                <th scope="col">Strand</th>
                <th width="12.5em" scope="col">Product</th>
                <th scope="col">Status</th>
                <th scope="col">Action</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(curr, index) in phageAnnotations" :key="index">
                <!-- ID -->
                <td>{{ curr.id }}</td>
                <!-- Left -->
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else>{{ curr.left }}</td>
                <!-- Right -->
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else>{{ curr.right }}</td>
                <!-- Strand -->
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else-if="curr.strand === '-'">Complementary</td>
                <td v-else>Direct</td>
                <!-- Product -->
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td
                  v-else-if="curr.function === 'None selected'"
                  style="font-size: 1.2em"
                >
                  {{ curr.function }}
                </td>
                <td
                  v-else-if="curr.function.length < 17"
                  style="font-size: 1.2em"
                >
                  {{ curr.function.substring(1, curr.function.length) }}
                </td>
                <td v-else style="font-size: 1.2em">
                  {{ curr.function.substring(1, 14) }}...
                </td>
                <!-- Status -->
                <td
                  v-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                ></td>
                <td v-else-if="curr.status === 'tRNA'" style="color: #7570b3">
                  {{ getStatus(index) }}
                </td>
                <td v-else-if="curr.status === 'Pass'" style="color: #1b9e77">
                  {{ getStatus(index) }}
                </td>
                <td v-else-if="curr.status === 'Fail'" style="color: #d95f02">
                  {{ getStatus(index) }}
                </td>
                <!-- Action -->
                <td
                  v-if="
                    (curr.function === '@DELETED' || curr.status === 'trnaDELETED') && !viewOnly
                  "
                >
                  <button
                    class="btn btn-outline-dark btn-sm"
                    style="width: 6.25em"
                    @click="reinstate(index)"
                  >
                    <strong>Reinstate</strong>
                  </button>
                </td>
                <td
                  v-else-if="
                    curr.function === '@DELETED' || curr.status === 'trnaDELETED'
                  "
                >
                </td>
                <td v-else-if="curr.status === 'tRNA' && !viewOnly">
                  <button
                    class="btn btn-outline-dark btn-sm"
                    style="width: 6.25em"
                    @click="deleteTRNA(index)"
                  >
                    <strong>Delete</strong>
                  </button>
                </td>
                <td v-else-if="curr.status === 'tRNA'">
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
                      style="width: 6.25em"
                      v-if="curr.function[0] === '@' && !viewOnly"
                    >
                      <strong>Annotate</strong>
                    </button>
                    <button
                      class="btn btn-outline-dark btn-sm"
                      style="width: 6.25em"
                      v-else-if="curr.function[0] === '@'"
                    >
                      <strong>View</strong>
                    </button>
                    <button
                      class="btn btn-dark btn-sm"
                      style="width: 6.25em"
                      v-else-if="curr.function === 'None selected' && !viewOnly"
                    >
                      <strong>Annotate</strong>
                    </button>
                    <button
                      class="btn btn-dark btn-sm"
                      style="width: 6.25em"
                      v-else-if="curr.function === 'None selected'"
                    >
                      <strong>View</strong>
                    </button>
                  </router-link>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
      <div class="alert alert-secondary">
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Blast', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'GeneMap',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav" id="next-top">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
    <b-modal
      class="text-size"
      v-model="showAddCDS"
      ref="addCDSModal"
      id="addCDS-modal"
      hide-footer
    >
      <template #modal-title>
        <div class="text-size">Add CDS</div>
      </template>
      <b-form @submit="onSubmitAdd" align="left">
        <b-form-group label="Left:" label-size="lg" label-for="add-left-input">
          <b-form-input
            class="form-input"
            id="add-left-input"
            type="number"
            v-model="addCDS.left"
            required
            placeholder="Enter left position"
          ></b-form-input>
        </b-form-group>
        <b-form-group label="Right:" label-size="lg" label-for="add-right-input">
          <b-form-input
            class="form-input"
            id="add-right-input"
            type="number"
            v-model="addCDS.right"
            required
            placeholder="Enter right position"
          ></b-form-input>
        </b-form-group>
        <b-form-group label="Strand:" label-size="lg">
          <b-form-select
            class="form-input"
            v-model="addCDS.strand"
            required
            :options="strandOptions"
          ></b-form-select>
        </b-form-group>
        <b-form-group
          label="Force Add:"
          label-size="lg"
          label-for="add-force-input"
        >
          <b-form-checkbox
            id="add-force-input"
            type="checkbox"
            v-model="addCDS.force"
          >
            Do not check this unless you are sure you want to add the CDS. This
            could remove another CDS if a duplicate ID is given or create a CDS
            that does not represent an actual gene.
          </b-form-checkbox>
        </b-form-group>
        <hr />
        <b-button type="submit" class="mt-3" block style="margin-top: 0em">
          <strong>Submit</strong>
        </b-button>
      </b-form>
    </b-modal>
    <b-toast id="annotations-status" variant="primary" no-auto-hide>
      <template #toast-title>
        <strong class="text-size"> {{statusTitle}} </strong>
      </template>
      <div class="text-size">{{ statusMessage }}</div>
    </b-toast>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";
import { LoaderPlugin } from 'vue-google-login';
import Vue from 'vue';

export default {
  name: "Annotations",
  components: {
    Loading,
    Navbar,
  },

  data() {
    return {
      viewOnly: false,
      addCDS: {
        id: "",
        left: "",
        right: "",
        strand: null,
        force: false,
        read: [],
      },
      strandOptions: [
        { value: null, text: "Please select a direction" },
        { value: "+", text: "+ (Direct)" },
        { value: "-", text: "- (Complementary)" },
      ],
      phageAnnotations: [],
      currCDS: {
        id: "",
        left: "",
        right: "",
        strand: "",
      },
      showAddCDS: false,
      pageLoading: true,
      blastLoading: false,
      completedGenes: 0,
      gap: 10,
      overlap: 10,
      oppositeGap: 50,
      short: 200,
      interval: null,
      waitMessage: "Your Blast results will be interpreted.",
      statusMessage: "",
      statusTitle: "",
    };
  },

  beforeCreate() {
    Vue.GoogleAuth.then(auth2 => {
      if (!auth2.isSignedIn.get()) {
        this.$router.push('/');
      }
      axios
        .get(process.env.VUE_APP_BASE_URL + `/check_user/${auth2.currentUser.get().ft.Qt}/${this.$route.params.phageID}`)
        .then((response) => {
          if (response.data === "fail") {
            this.$router.push('/');
          }
          else if (response.data.view) {
            this.viewOnly = true
          }
        })
        .catch((error) => {
          console.error(error);
        });
    })
  },

  created() {
    this.getData();
  },

  destroyed() {
    this.stopChecking();
  },

  computed: {
    navUpload: function () {
      return !this.viewOnly;
    },

    navBlast: function () {
      return !this.viewOnly;
    },

    navAnnotations: function () {
      return true;
    },

    navGeneMap: function () {
      return true;
    },

    navSettings: function () {
      return true;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
    },
  },

  methods: {
    /**
     * Populates phageAnnotations with data.
     */
    getData() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/none`
        )
        .then((response) => {
          this.phageAnnotations = response.data.annotations;
          console.log(this.phageAnnotations[0]);
          console.log(this.phageAnnotations[1]);
          this.gap = response.data.gap;
          this.overlap = response.data.overlap;
          this.oppositeGap = response.data.opposite_gap;
          this.short = response.data.short;
          this.pageLoading = false;
          for (var i = 0; i < this.phageAnnotations.length; i += 1) {
            if (this.phageAnnotations[i].function !== 'None selected')
              this.completedGenes += 1;
          }
          this.parseBlast();
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Calls method on back-end that parses and stores all of the BLAST results.
     */
    parseBlast() {
      axios({
        method: 'get',
        url: process.env.VUE_APP_BASE_URL + `/annotations/${this.$route.params.phageID}/blast`,
      })
        .then((response) => {
          console.log(response.data);
          if (response.data === "empty") {
            this.blastLoading = true;
            this.checkIfParseComplete();
          }
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Stops the interval loop for checkIfFilesUploaded().
     */
    stopChecking() {
      clearInterval(this.interval);
    },

    /**
     * Checks to see if the fasta file has been uploaded.
     */
    checkIfParseComplete() {
      this.interval = setInterval(() => {
        axios
          .get(
            process.env.VUE_APP_BASE_URL +
              `/annotations/${this.$route.params.phageID}/check`
          )
          .then((response) => {
            console.log(response.data);
            if (response.data === 'success') {
              this.blastLoading = false;
              this.stopChecking();
              this.statusTitle = "FINISHED";
              this.statusMessage = "All of the BLAST results have finished being interpreted.";
              this.$bvToast.show('annotations-status');
            } else if (response.data === 'error') {
              this.blastLoading = false;
              this.stopChecking();
              this.statusTitle = "ERROR";
              this.statusMessage = `An unknown error occurred. 
                                    Try removing and reuploading the BLAST files. 
                                    If you ignore this error not all of your BLAST results or ORFs will be shown. 
                                    If this error continues, please contact us by visiting the 'about' page.`;
              this.$bvToast.show('annotations-status');
            } else if (response.data !== 'complete') {
              if (response.data === '0') {
                this.waitMessage = "Your Blast results are currently being interpreted. This may take several minutes."
              } else {
                this.waitMessage = "Your Blast results are number " + response.data + " in line to be interpreted."
              }
            }
          })
          .catch((error) => {
            console.error(error);
          });
      }, 5000);
    },

    /**
     * Changes the function of a deleted gene so that it is now visible.
     */
    reinstate(index) {
      if (this.phageAnnotations[index].function === '@DELETED') {
        this.completedGenes -= 1;
        this.phageAnnotations[index].function = 'None selected';
      } else {
        this.phageAnnotations[index].status = 'tRNA';
      }
      const payload = {
        id: this.phageAnnotations[index].id,
        left: this.phageAnnotations[index].left,
        right: this.phageAnnotations[index].right,
        strand: this.phageAnnotations[index].strand,
        function: this.phageAnnotations[index].function,
        status: this.phageAnnotations[index].status,
        frame: this.phageAnnotations[index].frame,
      };
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${this.phageAnnotations[index].id}`,
          payload
        )
        .then(() => {
          console.log(this.phageAnnotations[index].function);
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Temporarily deletes a tRNA gene.
     */
    deleteTRNA(index) {
      this.phageAnnotations[index].status = 'trnaDELETED';
      const payload = {
        id: this.phageAnnotations[index].id,
        left: this.phageAnnotations[index].left,
        right: this.phageAnnotations[index].right,
        strand: this.phageAnnotations[index].strand,
        function: this.phageAnnotations[index].function,
        status: 'trnaDELETED',
        frame: this.phageAnnotations[index].frame,
      };
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${this.phageAnnotations[index].id}`,
          payload
        )
        .then(() => {
          console.log(this.phageAnnotations[index].function);
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Creates the status for a given cds.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getStatus(index) {
      if (this.phageAnnotations[index].status === 'tRNA') {
        return 'tRNA';
      }
      var status = '';
      if (
        this.phageAnnotations[index].right - this.phageAnnotations[index].left <
        this.short
      ) {
        status += 'S | ';
      }
      if (index === 0) {
        status += this.getTailingOverlap(index);
      } else if (index === this.phageAnnotations.length - 1) {
        status += this.getLeadingOverlap(index);
      } else {
        status += this.getLeadingOverlap(index);
        status += this.getTailingOverlap(index);
      }
      if (status === '') {
        return 'R';
      }
      return status.substring(0, status.length - 3);
    },

    /**
     * Determines if there is a large gap or overlap between the current gene and the next.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getTailingOverlap(index) {
      var status = '';
      var nextGene = 1;
      var overlap = -10000;
      while (overlap === -10000) {
        if (index + nextGene >= this.phageAnnotations.length) {
          return status;
        }
        if (this.phageAnnotations[index + nextGene].function !== '@DELETED') {
          overlap =
            this.phageAnnotations[index].right - this.phageAnnotations[index + nextGene].left;
        } else {
          nextGene += 1;
        }
      }
      if (
        this.phageAnnotations[index].strand === this.phageAnnotations[index + nextGene].strand
      ) {
        if (overlap < this.gap * -1) {
          status += 'LTG | ';
        } else if (overlap > this.overlap) {
          status += 'LTO | ';
        }
      } else if (overlap > this.oppositeGap * -1) {
        status += 'STG | ';
      }
      return status;
    },

    /**
     * Determines if there is a large gap or overlap between the current gene and the previous.
     * @param {number} index the index of a cds.
     * @return {string} the status of the cds.
     */
    getLeadingOverlap(index) {
      var status = '';
      var nextGene = 1;
      var leadingOverlap = -10000;
      while (leadingOverlap === -10000) {
        if (index - nextGene < 0) {
          return status;
        }
        if (this.phageAnnotations[index - nextGene].function !== '@DELETED') {
          leadingOverlap =
            this.phageAnnotations[index - nextGene].right - this.phageAnnotations[index].left;
        } else {
          nextGene += 1;
        }
      }
      if (
        this.phageAnnotations[index].strand === this.phageAnnotations[index - nextGene].strand
      ) {
        if (leadingOverlap < this.gap * -1) {
          status += 'LLG | ';
        } else if (leadingOverlap > this.overlap) {
          status += 'LLO | ';
        }
      } else if (leadingOverlap > this.oppositeGap * -1) {
        status += 'SLG | ';
      }
      return status;
    },

    /**
     * @param {event} evt the submit event on the add cds modal.
     * Adds a new CDS or notifies the user of a failed add.
     */
    onSubmitAdd(evt) {
      evt.preventDefault();
      this.$refs.addCDSModal.hide();
      let read = true;
      const payload = {
        id: 'added',
        left: this.addCDS.left,
        right: this.addCDS.right,
        strand: this.addCDS.strand,
        force: this.addCDS.force,
        read, // property shorthand
      };
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/none`,
          payload
        )
        .then((response) => {
          console.log(response.data.message);
          if (response.data.message === 'ID already exists.') {
            this.statusTitle = "ADD FAILED";
            this.statusMessage = `The CDS already exists. Try again with a different ORF. 
                                  To ignore this warning and add the CDS, check the 'Force Add' box.`;
            this.$bvToast.show('annotations-status');
          } else if (response.data.message === 'Not orf.') {
            this.statusTitle = "ADD FAILED";
            this.statusMessage = `The inputted left and right locations do not represent an ORF. 
                                  To ignore this warning and add the CDS, check the 'Force Add' box.`;
            this.$bvToast.show('annotations-status');
          } else {
            window.location.reload();
          }
        })
        .catch((error) => {
          console.error(error);
        });
    },
  },
};
</script>

<style scoped>
.wrapper {
  margin: 0;
}

h1 {
  margin-top: 0.7em;
}

.nav-btns-wrapper {
  text-align: center;
}

td,
th {
  word-wrap: break-word;
  width: 9.5em;
  font-size: 1.4em;
}

.table-responsive thead th {
  position: sticky;
  top: 0;
  background: white;
  color: black;
  border: darkgray;
  font-size: 1.4em;
}

.btn-nav {
  margin: 0.25em;
}

.btn-action {
  margin: 0.25em;
}

.btn-dark {
  font-size: 15pt;
}

.btn-outline-dark {
  font-size: 15pt;
}

.alert-secondary {
  background-color: white;
  border-color: white;
  font-size: 1.4em;
  text-align: left;
}

.text-size {
  font-size: 1.2em;
}

.form-input {
  height: 2em; 
  font-size: 15pt;
}
</style>
