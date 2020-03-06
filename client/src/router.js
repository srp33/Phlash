import Vue from 'vue';
import Router from 'vue-router';
import Home from './components/Home.vue';
import Ping from './components/Ping.vue';
import Upload from './components/Upload.vue';
import DNAMaster from './components/DNAMaster.vue';
import Upload_GM from './components/Upload_GM.vue';
import Annotation from './components/Annotation.vue';
import Failed from './components/Failed.vue';
import More from './components/More.vue';
import Upload_Test from './components/Upload_Test.vue';
import Blast from './components/Blast.vue';
import GenBank from './components/GenBank.vue';


Vue.use(Router);

export default new Router({
  mode: 'history',
  base: process.env.BASE_URL,
  routes: [
    {
      path: '/',
      name: 'Home',
      component: Home,
    },
    {
      path: '/ping',
      name: 'Ping',
      component: Ping,
    },
    {
      path: '/upload',
      name: 'Upload',
      component: Upload,
    },
    {
      path: '/dnamaster',
      name: 'DNAMaster',
      component: DNAMaster,
    },
    {
      path: '/api/upload_genemark',
      name: 'Upload_GM',
      component: Upload_GM,
    },
    {
      path: '/annotate_data',
      name: 'Annotation',
      component: Annotation,
    },
    {
      path: '/annotate_data/failed/:id',
      name: 'Failed',
      component: Failed,
    },
    {
      path: '/annotate_data/more/:id',
      name: 'More',
      component: More,
    },
    {
      path: '/test',
      name: 'Test',
      component: Upload_Test,
    },
    {
      path: '/blast/:id',
      name: 'Blast',
      component: Blast,
    },
    {
      path: '/genbank',
      name: 'GenBank',
      component: GenBank,
    },
  ],
});
