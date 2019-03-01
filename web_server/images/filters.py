from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin
from .tables import ImageTable
from .models import Image


class FilteredImageListView(SingleTableMixin, FilterView):
    table_class = ImageTable
    model = Image
    template_name = 'template.html'

    filterset_class = PersonFilter