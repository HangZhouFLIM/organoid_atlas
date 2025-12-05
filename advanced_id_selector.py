import os
import pandas as pd
import numpy as np
from qtpy.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QMessageBox,
    QFileDialog,
    QListWidget,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QWidget,
    QTextEdit,
)
from qtpy.QtCore import Qt


class AdvancedIDSelector(QDialog):
    def __init__(self, file_paths, parent=None):
        super().__init__(parent)
        self.file_paths = file_paths
        self.file_data = {}
        self.current_file = None
        self.setWindowTitle("Advanced ID Selector")
        self.setMinimumSize(1000, 700)
        self.setup_ui()
        self.load_files()

    def setup_ui(self):
        main_layout = QHBoxLayout()
        splitter = QSplitter(Qt.Horizontal)

        # File list widget
        self.file_list = QListWidget()
        self.file_list.itemClicked.connect(self.file_selected)
        splitter.addWidget(self.file_list)

        # Right side panel
        right_panel = QVBoxLayout()

        # File info
        self.file_label = QLabel("No file selected")
        right_panel.addWidget(self.file_label)

        # Tab widget for different functions
        self.tab_widget = QTabWidget()

        # Tab 1: Select IDs
        self.select_tab = QWidget()
        self.setup_select_tab()
        self.tab_widget.addTab(self.select_tab, "Select IDs")

        # Tab 2: Exclude IDs
        self.exclude_tab = QWidget()
        self.setup_exclude_tab()
        self.tab_widget.addTab(self.exclude_tab, "Exclude IDs")

        # Tab 3: Preview Data
        self.preview_tab = QWidget()
        self.setup_preview_tab()
        self.tab_widget.addTab(self.preview_tab, "Preview Data")

        right_panel.addWidget(self.tab_widget)

        # Export buttons
        buttons_layout = QHBoxLayout()

        self.export_button = QPushButton("Export All Files")
        self.export_button.clicked.connect(self.export_all_data)
        buttons_layout.addWidget(self.export_button)

        self.combined_export_button = QPushButton("Export Combined Data")
        self.combined_export_button.clicked.connect(self.export_combined_data)
        buttons_layout.addWidget(self.combined_export_button)

        # Add new button for NaN/Inf stats
        self.stats_button = QPushButton("Show Data Quality Report")
        self.stats_button.clicked.connect(self.show_data_quality_report)
        buttons_layout.addWidget(self.stats_button)

        right_panel.addLayout(buttons_layout)

        right_widget = QWidget()
        right_widget.setLayout(right_panel)
        splitter.addWidget(right_widget)

        main_layout.addWidget(splitter)
        self.setLayout(main_layout)

    def setup_select_tab(self):
        layout = QVBoxLayout()

        self.select_label = QLabel("Enter IDs to SELECT (comma separated):")
        layout.addWidget(self.select_label)

        self.select_input = QLineEdit()
        layout.addWidget(self.select_input)

        self.select_button = QPushButton("Set Selected IDs")
        self.select_button.clicked.connect(self.process_selected_ids)
        layout.addWidget(self.select_button)

        self.current_selected_label = QLabel("Currently selected IDs: None")
        layout.addWidget(self.current_selected_label)

        self.select_tab.setLayout(layout)

    def setup_exclude_tab(self):
        layout = QVBoxLayout()

        self.exclude_label = QLabel("Enter IDs to EXCLUDE from unselected (comma separated):")
        layout.addWidget(self.exclude_label)

        self.exclude_input = QLineEdit()
        layout.addWidget(self.exclude_input)

        self.exclude_button = QPushButton("Set Excluded IDs")
        self.exclude_button.clicked.connect(self.process_excluded_ids)
        layout.addWidget(self.exclude_button)

        self.current_excluded_label = QLabel("Currently excluded IDs: None")
        layout.addWidget(self.current_excluded_label)

        self.exclude_tab.setLayout(layout)

    def setup_preview_tab(self):
        layout = QVBoxLayout()

        self.preview_table = QTableWidget()
        self.preview_table.setColumnCount(5)
        self.preview_table.setHorizontalHeaderLabels(
            ["ID", "MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
        )
        self.preview_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        layout.addWidget(self.preview_table)
        self.preview_tab.setLayout(layout)

    def load_files(self):
        for file_path in self.file_paths:
            try:
                df = pd.read_csv(file_path)
                if not df.empty:
                    filename = os.path.basename(file_path)
                    self.file_data[filename] = {
                        "path": file_path,
                        "df": df,
                        "selected_ids": [],
                        "excluded_ids": [],
                    }
                    self.file_list.addItem(filename)
            except Exception as e:  # noqa: BLE001
                print(f"Error loading {file_path}: {str(e)}")

        if self.file_list.count() > 0:
            self.file_list.setCurrentRow(0)
            self.file_selected(self.file_list.currentItem())

    def file_selected(self, item):
        filename = item.text()
        self.current_file = filename
        file_info = self.file_data[filename]
        df = file_info["df"]

        # Update file info
        self.file_label.setText(
            f"Selected: {filename}\n"
            f"Total IDs: {len(df)}\n"
            f"Selected IDs: {len(file_info['selected_ids'])}\n"
            f"Excluded IDs: {len(file_info['excluded_ids'])}"
        )

        # Update select tab
        self.select_input.setText(", ".join(map(str, file_info["selected_ids"])))
        self.current_selected_label.setText(
            f"Currently selected IDs: {', '.join(map(str, file_info['selected_ids'])) or 'None'}"
        )

        # Update exclude tab
        self.exclude_input.setText(", ".join(map(str, file_info["excluded_ids"])))
        self.current_excluded_label.setText(
            f"Currently excluded IDs: {', '.join(map(str, file_info['excluded_ids'])) or 'None'}"
        )

        # Update preview
        self.update_preview_table()

    def update_preview_table(self):
        if not self.current_file:
            return

        file_info = self.file_data[self.current_file]
        df = file_info["df"]

        # Find available columns (case insensitive)
        target_columns = ["MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
        available_columns = []

        for col in df.columns:
            for target in target_columns:
                if target.lower() in col.lower():
                    available_columns.append(col)
                    break

        if not available_columns:
            self.preview_table.setRowCount(0)
            return

        # Get all IDs
        all_ids = df.iloc[:, 0].unique()

        # Filter based on selection and exclusion
        selected_ids = set(file_info["selected_ids"])
        excluded_ids = set(file_info["excluded_ids"])

        display_ids = []
        for id_val in all_ids:
            if id_val in selected_ids:
                display_ids.append((id_val, "Selected"))
            elif id_val in excluded_ids:
                display_ids.append((id_val, "Excluded"))
            else:
                display_ids.append((id_val, "Unselected"))

        # Setup table
        self.preview_table.setRowCount(len(display_ids))
        self.preview_table.setColumnCount(len(available_columns) + 2)
        headers = ["ID", "Status"] + available_columns
        self.preview_table.setHorizontalHeaderLabels(headers)

        for row, (id_val, status) in enumerate(display_ids):
            # Add ID and status
            self.preview_table.setItem(row, 0, QTableWidgetItem(str(id_val)))
            self.preview_table.setItem(row, 1, QTableWidgetItem(status))

            # Add data columns
            row_data = df[df.iloc[:, 0] == id_val].iloc[0]
            for col_idx, col_name in enumerate(available_columns, start=2):
                value = str(row_data[col_name]) if col_name in row_data else "N/A"
                self.preview_table.setItem(row, col_idx, QTableWidgetItem(value))

    def process_selected_ids(self):
        if not self.current_file:
            return

        input_text = self.select_input.text()
        file_info = self.file_data[self.current_file]
        df = file_info["df"]

        try:
            ids = [int(id_str.strip()) for id_str in input_text.split(",") if id_str.strip()]
            valid_ids = [id_val for id_val in ids if id_val in df.iloc[:, 0].values]

            if ids and not valid_ids:
                QMessageBox.warning(self, "Warning", "No valid IDs found in the dataset!")
                return

            file_info["selected_ids"] = valid_ids
            self.current_selected_label.setText(
                f"Currently selected IDs: {', '.join(map(str, valid_ids)) or 'None'}"
            )

            msg = f"Selected {len(valid_ids)} IDs for {self.current_file}"
            if len(valid_ids) != len(ids):
                msg += f" ({len(ids) - len(valid_ids)} invalid IDs ignored)"
            QMessageBox.information(self, "Success", msg)

            self.update_preview_table()
            self.file_selected(self.file_list.currentItem())
        except ValueError:
            QMessageBox.warning(self, "Error", "Please enter valid integer IDs separated by commas")

    def process_excluded_ids(self):
        if not self.current_file:
            return

        input_text = self.exclude_input.text()
        file_info = self.file_data[self.current_file]
        df = file_info["df"]

        try:
            ids = [int(id_str.strip()) for id_str in input_text.split(",") if id_str.strip()]
            valid_ids = [id_val for id_val in ids if id_val in df.iloc[:, 0].values]

            # Check if any IDs are already selected
            selected_ids = set(file_info["selected_ids"])
            conflicting_ids = [id_val for id_val in valid_ids if id_val in selected_ids]

            if conflicting_ids:
                QMessageBox.warning(
                    self,
                    "Conflict",
                    "These IDs are already selected and cannot be excluded: "
                    f"{', '.join(map(str, conflicting_ids))}",
                )
                valid_ids = [id_val for id_val in valid_ids if id_val not in selected_ids]

            if ids and not valid_ids:
                QMessageBox.warning(self, "Warning", "No valid IDs found to exclude!")
                return

            file_info["excluded_ids"] = valid_ids
            self.current_excluded_label.setText(
                f"Currently excluded IDs: {', '.join(map(str, valid_ids)) or 'None'}"
            )

            msg = f"Excluded {len(valid_ids)} IDs for {self.current_file}"
            if len(valid_ids) != len(ids):
                msg += f" ({len(ids) - len(valid_ids)} invalid IDs ignored)"
            QMessageBox.information(self, "Success", msg)

            self.update_preview_table()
            self.file_selected(self.file_list.currentItem())
        except ValueError:
            QMessageBox.warning(self, "Error", "Please enter valid integer IDs separated by commas")

    def export_all_data(self):
        if not self.file_data:
            QMessageBox.warning(self, "Warning", "No files loaded!")
            return

        # Get export directory
        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", "")

        if not export_dir:
            return

        exported_files = []

        for filename, file_info in self.file_data.items():
            df = file_info["df"]
            selected_ids = file_info["selected_ids"]
            excluded_ids = file_info["excluded_ids"]
            original_path = file_info["path"]
            base_name = os.path.splitext(os.path.basename(original_path))[0]

            # Get the required columns (case insensitive)
            target_columns = ["MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
            available_columns = []

            for col in df.columns:
                for target in target_columns:
                    if target.lower() in col.lower():
                        available_columns.append(col)
                        break

            if not available_columns:
                print(f"Skipping {filename} - required columns not found")
                continue

            # Export selected data
            if selected_ids:
                selected_data = df[df.iloc[:, 0].isin(selected_ids)]
                selected_output = os.path.join(export_dir, f"{base_name}_selected.csv")
                selected_data[available_columns].to_csv(selected_output, index=False)
                exported_files.append(selected_output)

            # Export unselected data (excluding the excluded IDs)
            unselected_mask = ~df.iloc[:, 0].isin(selected_ids)
            if excluded_ids:
                unselected_mask &= ~df.iloc[:, 0].isin(excluded_ids)

            unselected_data = df[unselected_mask]
            unselected_output = os.path.join(export_dir, f"{base_name}_unselected.csv")
            unselected_data[available_columns].to_csv(unselected_output, index=False)
            exported_files.append(unselected_output)

            # Export excluded data if any
            if excluded_ids:
                excluded_data = df[df.iloc[:, 0].isin(excluded_ids)]
                excluded_output = os.path.join(export_dir, f"{base_name}_excluded.csv")
                excluded_data[available_columns].to_csv(excluded_output, index=False)
                exported_files.append(excluded_output)

        if exported_files:
            QMessageBox.information(
                self,
                "Success",
                f"Exported {len(exported_files)} files to:\n{export_dir}",
            )
        else:
            QMessageBox.warning(self, "Warning", "No data was exported!")

    def export_combined_data(self):
        if not self.file_data:
            QMessageBox.warning(self, "Warning", "No files loaded!")
            return

        # Get export directory
        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", "")

        if not export_dir:
            return

        # Initialize combined DataFrames
        combined_selected = pd.DataFrame()
        combined_unselected = pd.DataFrame()

        for filename, file_info in self.file_data.items():
            df = file_info["df"]
            selected_ids = file_info["selected_ids"]
            excluded_ids = file_info["excluded_ids"]

            # Get available columns (case insensitive)
            target_columns = ["MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
            available_columns = [
                col for col in df.columns if any(c.lower() in col.lower() for c in target_columns)
            ]

            if not available_columns:
                continue

            # Add filename as a column
            df["SourceFile"] = filename

            # Process selected data
            if selected_ids:
                selected_data = df[df.iloc[:, 0].isin(selected_ids)]
                selected_data["Group"] = "Selected"
                combined_selected = pd.concat([combined_selected, selected_data])

            # Process unselected data (excluding excluded IDs)
            unselected_mask = ~df.iloc[:, 0].isin(selected_ids)
            if excluded_ids:
                unselected_mask &= ~df.iloc[:, 0].isin(excluded_ids)

            unselected_data = df[unselected_mask]
            if not unselected_data.empty:
                unselected_data["Group"] = "Unselected"
                combined_unselected = pd.concat([combined_unselected, unselected_data])

        # Export combined files if they contain data
        exported_files = []

        if not combined_selected.empty:
            selected_output = os.path.join(export_dir, "combined_selected.csv")
            combined_selected.to_csv(selected_output, index=False)
            exported_files.append(selected_output)

        if not combined_unselected.empty:
            unselected_output = os.path.join(export_dir, "combined_unselected.csv")
            combined_unselected.to_csv(unselected_output, index=False)
            exported_files.append(unselected_output)

        if exported_files:
            QMessageBox.information(
                self,
                "Success",
                f"Exported combined files to:\n{export_dir}",
            )
        else:
            QMessageBox.warning(self, "Warning", "No data was exported for combined files!")

    def show_data_quality_report(self):
        """Show comprehensive data quality report with NaN/Inf stats and cleaned data."""
        if not self.file_data:
            QMessageBox.warning(self, "Warning", "No files loaded!")
            return

        # Create a new dialog for the report
        report_dialog = QDialog(self)
        report_dialog.setWindowTitle("Data Quality Report")
        report_dialog.setMinimumSize(800, 600)

        layout = QVBoxLayout()

        # Add a text area for the report
        self.report_text = QTextEdit()
        self.report_text.setReadOnly(True)
        layout.addWidget(self.report_text)

        # Add a table for cleaned data preview
        self.cleaned_data_table = QTableWidget()
        self.cleaned_data_table.setColumnCount(5)
        self.cleaned_data_table.setHorizontalHeaderLabels(
            ["ID", "MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
        )
        self.cleaned_data_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        layout.addWidget(self.cleaned_data_table)

        # Add close button
        close_button = QPushButton("Close")
        close_button.clicked.connect(report_dialog.close)
        layout.addWidget(close_button)

        report_dialog.setLayout(layout)

        # Generate and display the report
        self._generate_data_quality_report()

        report_dialog.exec_()

    def _generate_data_quality_report(self):
        """Generate the data quality report with stats and cleaned data."""
        # Initialize combined DataFrames
        combined_selected = pd.DataFrame()
        combined_unselected = pd.DataFrame()

        for filename, file_info in self.file_data.items():
            df = file_info["df"]
            selected_ids = file_info["selected_ids"]
            excluded_ids = file_info["excluded_ids"]

            # Get available columns (case insensitive)
            target_columns = ["MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
            available_columns = [
                col for col in df.columns if any(c.lower() in col.lower() for c in target_columns)
            ]

            if not available_columns:
                continue

            # Process selected data
            if selected_ids:
                selected_data = df[df.iloc[:, 0].isin(selected_ids)]
                combined_selected = pd.concat([combined_selected, selected_data])

            # Process unselected data (excluding excluded IDs)
            unselected_mask = ~df.iloc[:, 0].isin(selected_ids)
            if excluded_ids:
                unselected_mask &= ~df.iloc[:, 0].isin(excluded_ids)

            unselected_data = df[unselected_mask]
            if not unselected_data.empty:
                combined_unselected = pd.concat([combined_unselected, unselected_data])

        # Prepare report text
        report_text = "DATA QUALITY REPORT\n\n"

        if not combined_selected.empty:
            report_text += "=== SELECTED DATA ===\n"
            report_text += self._get_data_stats(combined_selected, "Selected")
            cleaned_selected = self._get_cleaned_data(combined_selected)
            self._display_cleaned_data(cleaned_selected, "Selected")

        if not combined_unselected.empty:
            report_text += "\n=== UNSELECTED DATA ===\n"
            report_text += self._get_data_stats(combined_unselected, "Unselected")
            cleaned_unselected = self._get_cleaned_data(combined_unselected)
            self._display_cleaned_data(cleaned_unselected, "Unselected")

        if combined_selected.empty and combined_unselected.empty:
            report_text += "No data available to generate report."

        self.report_text.setPlainText(report_text)

    def _get_data_stats(self, df, group_name):
        """Generate statistics text for a DataFrame."""
        target_columns = ["MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
        available_columns = [
            col for col in df.columns if any(c.lower() in col.lower() for c in target_columns)
        ]

        if not available_columns:
            return f"No relevant columns found in {group_name} data.\n"

        stats = f"{group_name} Data Summary:\n"

        stats += f"Total IDs collected: {len(df)}\n\n"

        stats += "Invalid Values Count:\n"

        for col in available_columns:
            # Count NaN values
            nan_count = df[col].isna().sum()

            # Count infinite values
            try:
                inf_count = np.isinf(df[col]).sum()
            except TypeError:
                inf_count = 0  # non-numeric columns

            stats += f"{col}:\n"
            stats += f"  NaN values: {nan_count}\n"
            stats += f"  Inf values: {inf_count}\n"
            stats += f"  Valid values: {len(df) - nan_count - inf_count}\n"

        return stats

    def _get_cleaned_data(self, df):
        """Return a cleaned DataFrame with NaN/Inf values removed."""
        target_columns = ["MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
        available_columns = [
            col for col in df.columns if any(c.lower() in col.lower() for c in target_columns)
        ]

        if not available_columns:
            return pd.DataFrame()

        # Make a copy to avoid modifying original
        cleaned_df = df.copy()

        # Replace inf with NaN first
        for col in available_columns:
            try:
                cleaned_df[col] = cleaned_df[col].replace([np.inf, -np.inf], np.nan)
            except TypeError:
                pass  # skip non-numeric columns

        # Drop rows with any NaN values in our target columns
        cleaned_df = cleaned_df.dropna(subset=available_columns)

        return cleaned_df

    def _display_cleaned_data(self, cleaned_df, group_name):
        """Display cleaned data in the table."""
        if cleaned_df.empty:
            return

        target_columns = ["MeanInt", "MeanTauPhase", "MeanTauModulation", "Circularity"]
        available_columns = [
            col for col in cleaned_df.columns if any(c.lower() in col.lower() for c in target_columns)
        ]

        # Get ID column (assuming it's the first column)
        id_col = cleaned_df.columns[0]

        # Setup table
        self.cleaned_data_table.setRowCount(len(cleaned_df))
        self.cleaned_data_table.setColumnCount(len(available_columns) + 1)
        headers = [id_col] + available_columns
        self.cleaned_data_table.setHorizontalHeaderLabels(headers)

        # Populate table
        for row_idx, (_, row_data) in enumerate(cleaned_df.iterrows()):
            # Add ID
            self.cleaned_data_table.setItem(row_idx, 0, QTableWidgetItem(str(row_data[id_col])))

            # Add data columns
            for col_idx, col_name in enumerate(available_columns, start=1):
                value = str(row_data[col_name]) if col_name in row_data else "N/A"
                self.cleaned_data_table.setItem(row_idx, col_idx, QTableWidgetItem(value))


def run_advanced_selector():
    # Open file dialog to select multiple CSVs
    file_paths, _ = QFileDialog.getOpenFileNames(None, "Select CSV Files", "", "CSV Files (*.csv)")

    if not file_paths:
        return

    dialog = AdvancedIDSelector(file_paths)
    dialog.exec_()


# Run the function
if __name__ == "__main__":
    run_advanced_selector()
