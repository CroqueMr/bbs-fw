using BBSFW.ViewModel.Base;
using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Input;

namespace BBSFW.ViewModel
{
        public class SystemViewModel : ObservableObject
        {

                private ConfigurationViewModel _configVm;
                public ConfigurationViewModel ConfigVm
                {
                        get { return _configVm; }
                }

                private ConnectionViewModel _connectionVm;
                public ConnectionViewModel ConnectionVm
                {
                        get { return _connectionVm; }
                }

                public static List<ValueItemViewModel<bool>> OperationModeOptions { get; } =
                        new List<ValueItemViewModel<bool>>
                        {
                                new ValueItemViewModel<bool>(false, "Legal"),
                                new ValueItemViewModel<bool>(true, "Offroad")
                        };

                private ValueItemViewModel<bool> _selectedOperationMode;
                public ValueItemViewModel<bool> SelectedOperationMode
                {
                        get { return _selectedOperationMode; }
                        set
                        {
                                if (_selectedOperationMode != value)
                                {
                                        _selectedOperationMode = value;
                                        OnPropertyChanged(nameof(SelectedOperationMode));
                                }
                        }
                }

                public ICommand SetOperationModeCommand
                {
                        get { return new DelegateCommand(OnSetOperationMode); }
                }

                public SystemViewModel(ConfigurationViewModel config, ConnectionViewModel connectionVm)
                {
                        _configVm = config;
                        _connectionVm = connectionVm;
                        SelectedOperationMode = OperationModeOptions[0];
                }

                private async void OnSetOperationMode()
                {
                        if (!_connectionVm.IsConnected)
                        {
                                MessageBox.Show("Not Connected!", "Error", MessageBoxButton.OK, MessageBoxImage.Error);
                                return;
                        }

                        var res = await _connectionVm.GetConnection().SetOperationMode(SelectedOperationMode.Value, TimeSpan.FromSeconds(3));
                        if (!res.Timeout)
                        {
                                if (res.Result)
                                {
                                        MessageBox.Show("Operation mode set!", "Success", MessageBoxButton.OK, MessageBoxImage.Information);
                                }
                                else
                                {
                                        MessageBox.Show("Failed to set operation mode, check log.", "Error", MessageBoxButton.OK, MessageBoxImage.Error);
                                }
                        }
                        else
                        {
                                MessageBox.Show("Failed to set operation mode, timeout occured.", "Error", MessageBoxButton.OK, MessageBoxImage.Error);
                        }
                }

        }
}
