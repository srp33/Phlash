import Vue from 'vue';
import Router from 'vue-router';
import Home from './views/Home.vue';
import Upload from './views/Upload.vue';
import DNAMaster from './views/DNAMaster.vue';
import Blast from './views/Blast.vue';
import Annotations from './views/Annotations.vue';
import Pass from './views/Pass.vue';
import Fail from './views/Fail.vue';
import More from './views/More.vue';

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
      path: '/upload/:phageID',
      name: 'Upload',
      component: Upload,
    },
    {
      path: '/dnamaster/:phageID',
      name: 'DNAMaster',
      component: DNAMaster,
    },
    {
      path: '/blast/:phageID',
      name: 'Blast',
      component: Blast
    },
    {
      path: '/annotations/:phageID',
      name: 'Annotations',
      component: Annotations,
    },
    {
      path: '/annotations/pass/:phageID/:cdsID',
      name: 'Pass',
      component: Pass,
    },
    {
      path: '/annotations/fail/:phageID/:cdsID',
      name: 'Fail',
      component: Fail,
    },
    {
      path: '/annotations/more/:phageID/:cdsID',
      name: 'More',
      component: More,
    },
  ],
});
