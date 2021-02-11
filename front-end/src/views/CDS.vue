<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
    />
    <div class="container">
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <div class="headers">
        <h1>ID: {{ $route.params.cdsID }}</h1>
        <h4>Left: {{ currentCDS.start }}</h4>
        <h4>Right: {{ currentCDS.stop }}</h4>
        <h4>Strand: {{ currentCDS.strand }}</h4>
        <h4>Frame: {{ this.frame }}</h4>
        <h4>Called by: {{ this.calledBy }}</h4>
      </div>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>
          Here you can make any necessary changes to this CDS by editing its start and stop sites and its function. 
          Before making any final decisions you should review and consider all of the information below including 
          the alternative open reading frames, the coding potential graphs, and the BLAST results. 
          Click 'Save' to save changes made to this CDS or you can click 'Delete' if you wish to remove this CDS. 
          Above you can see generic information for the current CDS including which auto-annotation programs called that gene.
        </p>
        <hr />
        <div style="float: right; width: 50%;">
          <p><strong>Notes:</strong></p>
          <b-form-textarea
            id="textarea"
            v-model="notes"
            placeholder="Original Glimmer call @bp 408 has strength..."
            rows="3"
            max-rows="10"
          ></b-form-textarea>
        </div>
        <p>
          <strong>Your selected open reading frame:</strong> {{ newStart }}-{{ newStop }} {{ newStrand }}
        </p>
        <p><strong>Your selected function:</strong> {{ displayFunction }}</p>
        <button
          type="button"
          class="btn btn-light btn-action"
          @click="editCDS"
        >
          <strong>Save</strong>
        </button>
        <button
          type="button"
          class="btn btn-light btn-action"
          @click="deleteCDS($route.params.cdsID)"
        >
          <strong>Delete</strong>
        </button>
        <hr />
        <div class="nav-btns-wrapper">
          <button
            type="button"
            class="btn btn-light btn-action"
            @click="navPrevCDS()"
          >
            <strong>&#129052; Prev CDS</strong>
          </button>
          <button
            type="button"
            class="btn btn-light btn-action"
            @click="navNextCDS()"
          >
            <strong>Next CDS &#129054;</strong>
          </button>
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Annotations', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-light btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
        </div>
      </div>
      <div class="coding-potential-table">
        <h4 style="text-align: center; margin: 20px">Alternative Open Reading Frames</h4>
        <div style="overflow: hidden;">
          <div class="table-responsive" style="float: left; width: 40%;">
            <strong> Direct Strand </strong>
            <table id="cp-table" class="table table-hover">
              <thead>
                <tr>
                  <th scope="col">Open Reading Frame</th>
                  <th scope="col">Action</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="(start, index) in dirStartOptions" :key="index">
                  <th>{{ start }}-{{ dirStopOptions[index] }}  +</th>
                  <td>
                    <button
                      v-if="start + dirStopOptions[index] !== currentCDS.start + currentCDS.stop"
                      class="btn btn-dark btn-sm"
                      @click="setORF(start, dirStopOptions[index], '+')"
                    >
                      <strong>Select</strong>
                    </button>
                  </td>
                </tr>
              </tbody>
            </table>
          </div>
          <div class="table-responsive" style="float: right; width: 40%">
            <strong> Complementary Strand </strong>
            <table id="cp-table" class="table table-hover">
              <thead>
                <tr>
                  <th scope="col">Open Reading Frame</th>
                  <th scope="col">Action</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="(start, index) in compStartOptions" :key="index">
                  <th>{{ start }}-{{ compStopOptions[index] }}  -</th>
                  <td>
                    <button
                      v-if="start + compStopOptions[index] !== currentCDS.start + currentCDS.stop"
                      class="btn btn-dark btn-sm"
                      @click="setORF(start, compStopOptions[index], '-')"
                    >
                      <strong>Select</strong>
                    </button>
                  </td>
                </tr>
              </tbody>
            </table>
          </div>
        </div>
      </div>
      <hr />
      <div class="coding-potential">
        <h4 style="text-align: center; margin: 40px; height: 100%">
          GeneMark's Coding Potential Per Frame
        </h4>
        <div class="alert alert-primary">
          <p><strong>Key</strong></p>
          <p>
            <strong class="red-text">Orange Dashed Line:</strong> The selected start and stop positions.<br />
            <strong class="blue-text">Blue Line:</strong> The coding potential in relation to base number.<br />
            <strong class="green-text">Green Line:</strong> A 0.75 coding potential reference line.<br />
            <strong class="grey-text">Purple Line:</strong> The previous gene's stop position and the next
            gene's start position.<br />
          </p>
        </div>
        <div class="coding-potential-graphs">
          <div v-if="dataExists">
            <Graphs
              :data1="data1"
              :data2="data2"
              :data3="data3"
              :data4="data4"
              :data5="data5"
              :data6="data6"
              :start="currentCDS.start"
              :stop="currentCDS.stop"
              :frame="frame"
              :prevStop="prevStop"
              :nextStart="nextStart"
            />
          </div>
          <div v-else>
            <p>Loading graphs...</p>
          </div>
        </div>
        <p class="graphs-caption">
          <em>Frames are counted 1-6 (direct 1-3 and complementary 4-6).</em>
        </p>
      </div>
      <hr />
      <h4 style="text-align: center; margin: 40px">BLAST Results for {{ newStart }}-{{ newStop }}</h4>
      <BlastResults
        :blastResults="allBlastResults[currentCDS.start.toString() + '-' + currentCDS.stop.toString() + '  ' + currentCDS.strand]"
        :allowSelect="true"
        @newFunction="setFunction"
      />
      <!-- <div class="blast-results">
        <h4 style="text-align: center; margin: 40px">BLAST Results</h4>
        <strong> Direct Strand </strong>
        <div v-if="dataExists" class="table-responsive2" id="accordion" style="float: center; width: 100%; margin: 1em;">
          <div class="card border-dark" v-for="key in dirBlastKeys" :key="key">
            <div class="card-body">
              <h4 class="mb-0">
                <button
                  class="btn btn-outline-dark btn-blast"
                  data-toggle="collapse"
                  aria-expanded="false"
                  v-bind:data-target="'#' + key"
                  v-bind:aria-controls="key"
                  v-if="key !== currentCDS.start.toString() + '-' + currentCDS.stop.toString() + '  ' + currentCDS.strand"
                >
                  <strong>{{ key }}</strong>
                </button>
                <button
                  class="btn btn-dark btn-blast"
                  data-toggle="collapse"
                  aria-expanded="false"
                  v-bind:data-target="'#' + key"
                  v-bind:aria-controls="key"
                  v-else
                >
                  <strong>{{ key }}</strong>
                </button>
              </h4>
            </div>
            <div v-bind:id="key" class="collapse" data-parent="#accordion">
              <div class="card-body">
                <BlastResults
                  v-if="key === currentCDS.start.toString() + '-' + currentCDS.stop.toString() + '  ' + currentCDS.strand"
                  :blastResults="dirBlastResults[key]"
                  :allowSelect="true"
                  @newFunction="setFunction"
                />
                <BlastResults
                  v-else
                  :blastResults="dirBlastResults[key]"
                  :allowSelect="false"
                  @newFunction="setFunction"
                />
              </div>
            </div>
          </div>
        </div>
        <strong> Complementary Strand </strong>
        <div v-if="dataExists" class="table-responsive2" id="accordion" style="float: center; width: 100%; margin: 1em;">
          <div class="card border-dark" v-for="key in compBlastKeys" :key="key">
            <div class="card-body">
              <h4 class="mb-0">
                <button
                  class="btn btn-outline-dark btn-blast"
                  data-toggle="collapse"
                  aria-expanded="false"
                  v-bind:data-target="'#' + key"
                  v-bind:aria-controls="key"
                  v-if="key !== currentCDS.start.toString() + '-' + currentCDS.stop.toString() + '  ' + currentCDS.strand"
                >
                  <strong>{{ key }}</strong>
                </button>
                <button
                  class="btn btn-dark btn-blast"
                  data-toggle="collapse"
                  aria-expanded="false"
                  v-bind:data-target="'#' + key"
                  v-bind:aria-controls="key"
                  v-else
                >
                  <strong>{{ key }}</strong>
                </button>
              </h4>
            </div>
            <div v-bind:id="key" class="collapse" data-parent="#accordion">
              <div class="card-body">
                <BlastResults
                  v-if="key === currentCDS.start.toString() + '-' + currentCDS.stop.toString() + '  ' + currentCDS.strand"
                  :blastResults="compBlastResults[key]"
                  :allowSelect="true"
                  @newFunction="setFunction"
                />
                <BlastResults
                  v-else
                  :blastResults="compBlastResults[key]"
                  :allowSelect="false"
                  @newFunction="setFunction"
                />
              </div>
            </div>
          </div>
        </div>
      </div> -->
      <hr />
      <div class="alert alert-primary">
        <div style="float: right; width: 50%;">
          <p><strong>Notes:</strong></p>
          <b-form-textarea
            id="textarea"
            v-model="notes"
            placeholder="Original Glimmer call @bp 408 has strength..."
            rows="3"
            max-rows="10"
          ></b-form-textarea>
        </div>
        <p>
          <strong>Your selected open reading frame:</strong> {{ newStart }}-{{ newStop }} {{ newStrand }}
        </p>
        <p><strong>Your selected function:</strong> {{ displayFunction }}</p>
        <button
          type="button"
          class="btn btn-light btn-action"
          @click="editCDS"
        >
          <strong>Save</strong>
        </button>
        <button
          type="button"
          class="btn btn-light btn-action"
          @click="deleteCDS($route.params.cdsID)"
        >
          <strong>Delete</strong>
        </button>
        <hr />
        <div class="nav-btns-wrapper">
          <button
            type="button"
            class="btn btn-light btn-action"
            @click="navPrevCDS()"
          >
            <strong>&#129052; Prev CDS</strong>
          </button>
          <button
            type="button"
            class="btn btn-light btn-action"
            @click="navNextCDS()"
          >
            <strong>Next CDS &#129054;</strong>
          </button>
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Annotations', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-light btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue";
import BlastResults from "../components/BlastResults.vue";
import Graphs from "../components/Graphs.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "CDS",
  components: {
    BlastResults,
    Graphs,
    Loading,
    Navbar,
  },

  data() {

    return {
      dirStartOptions: [],
      compStartOptions: [],
      dirStopOptions: [],
      compStopOptions: [],
      dirBlastResults: [],
      compBlastResults: [],
      allBlastResults: [],
      currentCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: "",
        frame: "",
      },
      updatedCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: "",
      },
      frame: null,
      newFunction: "None selected",
      newStart: null,
      newStop: null,
      newStrand: null,
      data1: [{ x: [], y: [] }],
      data2: [{ x: [], y: [] }],
      data3: [{ x: [], y: [] }],
      data4: [{ x: [], y: [] }],
      data5: [{ x: [], y: [] }],
      data6: [{ x: [], y: [] }],
      nextCDS: null,
      prevCDS: null,
      nextStart: null,
      prevStop: null,
      calledBy: "",
      glimmer: "",
      genemark: "",
      phanotate: "",
      notes: "",
      dataExists: false,
      pageLoading: true,
      showFunction: false,
      showStart: false,
      saved: true,
    };

  },

  created() {
    this.getData(this.$route.params.cdsID);
    window.scrollTo(0, 0);
  },

  computed: {

    dirBlastKeys: function () {
      return Object.keys(this.dirBlastResults);
    },

    compBlastKeys: function () {
      return Object.keys(this.compBlastResults);
    },

    navUpload: function () {
      return true;
    },

    navBlast: function () {
      return true;
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
     * Gets the GeneMark, Blast, and auto-annotation data for the given CDS.
     * @param {string} cdsID the ID of the CDS.
     */
    getData(cdsID) {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${cdsID}`
        )
        .then((response) => {
          if (response.data.message != "Finished") {
            this.$router.push(
              `/annotations/${this.$route.params.phageID}`
            );
          }
          this.currentCDS = response.data.cds;
          if (this.currentCDS.function != "DELETED") {
            this.newFunction = this.currentCDS.function;
          }
          this.dirBlastResults = response.data.dir_blast;
          this.compBlastResults = response.data.comp_blast;
          this.allBlastResults = response.data.all_blast;
          this.dirStartOptions = response.data.dir_start_options;
          this.dirStopOptions = response.data.dir_stop_options;
          this.compStartOptions = response.data.comp_start_options;
          this.compStopOptions = response.data.comp_stop_options;
          this.newStart = this.currentCDS.start;
          this.newStop = this.currentCDS.stop;
          this.newStrand = this.currentCDS.strand;
          this.prevStop = response.data.prev_stop;
          this.nextStart = response.data.next_start;
          this.data1 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_1,
            },
          ];
          this.data2 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_2,
            },
          ];
          this.data3 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_3,
            },
          ];
          this.data4 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_4,
            },
          ];
          this.data5 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_5,
            },
          ];
          this.data6 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_6,
            },
          ];
          this.dataExists = true;
          this.pageLoading = false;
          this.frame = this.currentCDS.frame;
          this.notes = this.currentCDS.notes;
          this.glimmer = response.data.glimmer;
          this.genemark = response.data.genemark;
          this.phanotate = response.data.phanotate;
          var called = false;
          var cds = this.currentCDS.start.toString() + '-' + this.currentCDS.stop.toString() + ' ' + this.currentCDS.strand;
          if (this.glimmer.indexOf(cds) > -1) {
            called = true;
            this.calledBy += "Glimmer, ";
          }
          if (this.genemark.indexOf(cds) > -1) {
            called = true;
            this.calledBy += "GeneMark, ";
          }
          if (this.phanotate.indexOf(cds) > -1) {
            called = true;
            this.calledBy += "Phanotate, ";
          }
          if (called) { this.calledBy = this.calledBy.substring(0, this.calledBy.length - 2); }
          else { this.calledBy = "None"}
          this.nextCDS = response.data.nextCDS;
          this.prevCDS = response.data.prevCDS;
        })
        .catch((error) => {
          console.error(error);
        });
    },

    navNextCDS() {
      var cont = true;
      if (!this.saved) {
        cont = confirm("Are you sure you want to continue without saving?");
      }
      if (cont == true) {
        if (this.nextCDS != "undefined") {
          this.$route.params.cdsID = this.nextCDS;
          this.$router.push(
            `/annotations/cds/${this.$route.params.phageID}/${this.$route.params.cdsID}`
          );
          window.location.reload();
        } else {
          this.$router.push(`/annotations/${this.$route.params.phageID}`);
        }
      }
    },

    navPrevCDS() {
      var cont = true;
      if (!this.saved) {
        cont = confirm("Are you sure you want to continue without saving?");
      }
      if (cont == true) {
        if (this.prevCDS != "undefined") {
          this.$route.params.cdsID = this.prevCDS;
          this.$router.push(
            `/annotations/cds/${this.$route.params.phageID}/${this.$route.params.cdsID}`
          );
          window.location.reload();
        } else {
          this.$router.push(`/annotations/${this.$route.params.phageID}`);
        }
      }
    },

    /**
     * Changes the CDS data to reflect the user's changes.
     */
    editCDS() {
      this.saved = true;
      this.updatedCDS = this.currentCDS;
      this.updatedCDS.start = this.newStart;
      this.updatedCDS.stop = this.newStop;
      this.updatedCDS.function = "@" + this.newFunction;
      const payload = {
        id: this.updatedCDS.id,
        start: this.updatedCDS.start,
        stop: this.updatedCDS.stop,
        strand: this.updatedCDS.strand,
        function: this.updatedCDS.function,
        notes: this.notes,
        status: this.currentCDS.status,
      };
      this.updateCDS(payload, this.updatedCDS.id);
    },
    
    /**
     * Updates the database with the new changes on the current cds.
     * @param {dictionary} payload the new CDS data.
     * @param {string} cdsID the ID of the CDS to be updated.
     */
    updateCDS(payload, cdsID) {
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${cdsID}`,
          payload
        )
        .then(() => {
          console.log(response);
        })
        .catch((error) => {
          console.error(error);
        });
    },

    /**
     * Changes the function of the CDS to deleted.
     * @param {string} cdsID the ID of the CDS to be deleted.
     */
    deleteCDS(cdsID) {
      this.newFunction = "DELETED";
      this.editCDS();
    },

    /**
     * Updates the function.
     * @param {string} function the user selected function.
     */
    setFunction(funct) {
      this.saved = false;
      let match = funct.match(/(.*)##(.*)/);
      this.displayFunction = match[1];
      this.newFunction = funct;
    },

    /**
     * Updates the start and stop postions to what the user selected.
     * @param {number} start the user selected start.
     * @param {number} stop the user selected stop.
     */
    setORF(start, stop, strand) {
      this.saved = false;
      if (start != this.newStart || stop != this.newStop) {
        this.dataExists = false;
        this.newFunction = "";
        this.displayFunction = "";
        this.newStop = stop;
        this.currentCDS.stop = stop;
        this.newStart = start;
        this.newStrand = strand;
        this.currentCDS.start = start;
        this.currentCDS.strand = strand;
        this.calledBy = "";
        var called = false;
        var cds = this.currentCDS.start.toString() + '-' + this.currentCDS.stop.toString() + ' ' + this.currentCDS.strand;
        if (this.glimmer.indexOf(cds) > -1) {
          called = true;
          this.calledBy += "Glimmer, ";
        }
        if (this.genemark.indexOf(cds) > -1) {
          called = true;
          this.calledBy += "GeneMark, ";
        }
        if (this.phanotate.indexOf(cds) > -1) {
          called = true;
          this.calledBy += "Phanotate, ";
        }
        if (called) { this.calledBy = this.calledBy.substring(0, this.calledBy.length - 2); }
        else { this.calledBy = "None"}
        this.frame = ((start + 2) % 3) + 1;
        if (this.currentCDS.strand == "-") {
          this.frame = ((this.currentCDS.stop + 2) % 3) + 4;
        }
        else this.frame = ((this.currentCDS.start + 2) % 3) + 1;
        console.log(this.frame);
        this.$nextTick().then(() => {
          this.dataExists = true;
        });
      }
    },

  },
};
</script>

<style scoped>
/* ----- Title and Headers ----- */

.blast-results {
  margin: 1em;
  margin-bottom: 3em;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 10px;
}

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
  margin: 7px;
}

.red-text {
  color: #d95f02;
}

.blue-text {
  color: rgb(41, 41, 228);
}

.green-text {
  color: #1b9e77;
}

.grey-text {
  color: #7570b3;
}

/* ----- Coding Potential ----- */
.coding-potential-table {
  margin: 50px auto;
  width: 80%;
}

.table-responsive {
  max-height: 250px;
  overflow-y: auto;
  display: inline-block;
}

.table-responsive2 {
  max-height: 400px;
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
}

tbody {
  width: 100%;
}

.coding-potential-graphs {
  height: 650px;
}

.graphs-caption {
  margin-top: 15px;
  text-align: left;
  color: grey;
}

/* ----- Blast Results ----- */
.btn-blast {
  width: 100%;
}

/* Responsive Design */
@media only screen and (max-width: 1200px) {
  .coding-potential-graphs {
    height: 1300px;
  }
}
</style>
