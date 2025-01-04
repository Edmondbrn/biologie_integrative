from PyQt6.QtCore import QObject, pyqtSignal
from distances import Distances

class WorkerDistances(QObject):
    progress_signal = pyqtSignal(int)   # signal pour la progression
    finished = pyqtSignal()

    def __init__(self, distances_obj, df_ref, df_splicing, comparison_couples, output_dir, output_basename, n_cores):
        super().__init__()
        self.distances_obj : Distances = distances_obj
        self.df_ref = df_ref
        self.df_splicing = df_splicing
        self.comparison_couples = comparison_couples
        self.output_dir = output_dir
        self.output_basename = output_basename
        self.n_cores = n_cores

    def run(self):
        """Méthode appelée dans un thread séparé."""
        self.distances_obj.parallel_start_manual(
            df_ref=self.df_ref,
            df_splicing=self.df_splicing,
            comparison_couples=self.comparison_couples,
            output_dir=self.output_dir,
            output_basename=self.output_basename,
            n_cores=self.n_cores,
            progress_callback=self._emit_progress  # <--- on passe cette fonction
        )
        self.finished.emit()

    def _emit_progress(self, value):
        """Appelé pendant la boucle imap, émet le signal de progression."""
        self.progress_signal.emit(value)